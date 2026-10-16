//! Genomic neighbours: an Unknown family whose copies sit next to the same
//! classified family is most often another part of that family's element,
//! split off by RepeatModeler (an LTR and its internal region, the ends and
//! body of a DNA transposon). On the Drosophila benchmark, with the labels of
//! 56 LTR families hidden, every one of the 36 adjacent to one family in at
//! least 30% of their copies took the right class from it, 34 from their own
//! element's internal region.

use std::{
    collections::HashMap,
    fs,
    io::{BufRead, BufReader},
    path::Path,
};

use serde::Serialize;

use super::{
    decide::{is_non_te, Confidence, Disposition, DispositionKind},
    rmout::RmHit,
    taxonomy::top_level,
};
use crate::Result;

#[derive(Clone, Debug, Serialize)]
pub struct NeighbourParams {
    /// Largest gap between two hits that still counts as adjacent.
    pub max_gap: u64,
    /// Hits shorter than this are ignored (fragments, short simple repeats).
    pub min_hit: u64,
    /// Fewest copies a family needs, and fewest adjacent to one partner.
    pub min_copies: usize,
    /// Share of a family's copies adjacent to one partner for a call, and for
    /// a medium-confidence one.
    pub min_share: f64,
    pub medium_share: f64,
}

impl Default for NeighbourParams {
    fn default() -> Self {
        Self {
            max_gap: 50,
            min_hit: 50,
            min_copies: 5,
            min_share: 0.3,
            medium_share: 0.5,
        }
    }
}

/// Copies of one family next to one partner family, by the partner's place
/// relative to the copy (in the copy's orientation).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Partner {
    pub upstream: usize,
    pub downstream: usize,
    /// Partners on the opposite strand (consensus orientation is arbitrary,
    /// so these are parts of the same element too).
    pub opposite: usize,
}

impl Partner {
    fn total(&self) -> usize {
        self.upstream + self.downstream + self.opposite
    }
}

#[derive(Clone, Debug, Default)]
pub struct Adjacency {
    pub copies: usize,
    pub partners: HashMap<String, Partner>,
}

/// Adjacency of every family's copies to the hits of other families.
pub fn adjacency(by_family: &HashMap<String, Vec<RmHit>>, p: &NeighbourParams) -> HashMap<String, Adjacency> {
    let mut by_contig: HashMap<&str, Vec<&RmHit>> = HashMap::new();
    for h in by_family.values().flatten() {
        if h.end - h.start >= p.min_hit {
            by_contig.entry(&h.contig).or_default().push(h);
        }
    }
    let mut out: HashMap<String, Adjacency> = HashMap::new();
    for hits in by_contig.values_mut() {
        hits.sort_by_key(|h| (h.start, h.end));
        for (i, h) in hits.iter().enumerate() {
            let a = out.entry(h.family.clone()).or_default();
            a.copies += 1;
            let left = i.checked_sub(1).map(|j| hits[j]);
            let right = hits.get(i + 1).copied();
            for (n, on_left) in [(left, true), (right, false)] {
                let Some(n) = n else { continue };
                let gap = if on_left {
                    h.start.saturating_sub(n.end)
                } else {
                    n.start.saturating_sub(h.end)
                };
                if gap > p.max_gap || n.family == h.family {
                    continue;
                }
                let e = a.partners.entry(n.family.clone()).or_default();
                if n.minus != h.minus {
                    e.opposite += 1;
                } else if on_left != h.minus {
                    e.upstream += 1;
                } else {
                    e.downstream += 1;
                }
            }
        }
    }
    out
}

/// Labels of partner families with homology behind them, from
/// `rep classify` tables: `family_id` (name before `#`) -> label.
pub fn read_partner_labels(paths: &[impl AsRef<Path>]) -> Result<HashMap<String, String>> {
    let mut labels = HashMap::new();
    for path in paths {
        let mut lines = BufReader::new(fs::File::open(path.as_ref())?).lines();
        let Some(header) = lines.next().transpose()? else { continue };
        let ix = |name: &str| header.split('\t').position(|c| c == name);
        let (Some(id), Some(label), Some(evidence)) = (ix("family_id"), ix("label"), ix("evidence")) else {
            continue;
        };
        for line in lines {
            let line = line?;
            let f = line.split('\t').collect::<Vec<_>>();
            let (Some(i), Some(l), Some(e)) = (f.get(id), f.get(label), f.get(evidence)) else {
                continue;
            };
            let homology = ["protein:", "nucleotide:", "domains:"].iter().any(|k| e.contains(k));
            if homology && !l.starts_with("Unknown") && !is_non_te(l) {
                labels.insert(i.split('#').next().unwrap_or(i).to_owned(), (*l).to_owned());
            }
        }
    }
    Ok(labels)
}

/// Give an `Unknown` the class of the labelled family its copies sit next to.
pub fn apply_neighbours(
    mut d: Disposition,
    a: Option<&Adjacency>,
    labels: &HashMap<String, String>,
    p: &NeighbourParams,
) -> Disposition {
    let Some(a) = a else { return d };
    if d.label != "Unknown" || a.copies < p.min_copies {
        return d;
    }
    let Some((name, partner, label)) = a
        .partners
        .iter()
        .filter_map(|(n, x)| labels.get(n).map(|l| (n, x, l)))
        .max_by_key(|(n, x, _)| (x.total(), std::cmp::Reverse(n.as_str())))
    else {
        return d;
    };
    let share = partner.total() as f64 / a.copies as f64;
    if partner.total() < p.min_copies || share < p.min_share {
        return d;
    }
    d.evidence.push(format!(
        "adjacent_to:{name}:{label}(share={share:.2},upstream={},downstream={},opposite={},copies={})",
        partner.upstream, partner.downstream, partner.opposite, a.copies
    ));
    // an LTR on both sides of the same strand of its partner's copies
    let both_sides = partner.upstream.min(partner.downstream) * 4 >= partner.total();
    if top_level(label) == "LTR" && both_sides {
        d.evidence.push(format!("ltr_of:{name}"));
    }
    d.kind = DispositionKind::KnownTe;
    d.label = label.clone();
    d.confidence = if share >= p.medium_share {
        Confidence::Medium
    } else {
        Confidence::Low
    };
    d
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(contig: &str, start: u64, end: u64, minus: bool, family: &str) -> RmHit {
        RmHit {
            contig: contig.into(),
            start,
            end,
            minus,
            family: family.into(),
            class: "Unknown".into(),
            rep_start: 1,
            rep_end: end - start,
            divergence: 10.0,
            score: 500,
            id: 1,
        }
    }

    fn unknown() -> Disposition {
        Disposition {
            kind: DispositionKind::Unresolved,
            label: "Unknown".into(),
            confidence: Confidence::Low,
            evidence: vec![],
            warnings: vec![],
        }
    }

    #[test]
    fn an_ltr_flanking_its_internal_region_takes_its_class() {
        // ten elements LTR-int-LTR, five on each strand, and two lone LTRs
        let mut all = Vec::new();
        for i in 0..10u64 {
            let c = format!("c{i}");
            let minus = i % 2 == 1;
            all.push(hit(&c, 1000, 1300, minus, "ltr"));
            all.push(hit(&c, 1300, 6300, minus, "int"));
            all.push(hit(&c, 6300, 6600, minus, "ltr"));
        }
        all.push(hit("solo1", 100, 400, false, "ltr"));
        all.push(hit("solo2", 100, 400, true, "ltr"));
        let mut by_family: HashMap<String, Vec<RmHit>> = HashMap::new();
        for h in all {
            by_family.entry(h.family.clone()).or_default().push(h);
        }
        let p = NeighbourParams::default();
        let adj = adjacency(&by_family, &p);
        let ltr = &adj["ltr"];
        assert_eq!(ltr.copies, 22);
        assert_eq!(
            ltr.partners["int"],
            Partner {
                upstream: 10,
                downstream: 10,
                opposite: 0
            }
        );
        let labels = HashMap::from([("int".to_owned(), "LTR/Gypsy".to_owned())]);
        let d = apply_neighbours(unknown(), adj.get("ltr"), &labels, &p);
        assert_eq!((d.label.as_str(), d.confidence), ("LTR/Gypsy", Confidence::Medium));
        assert!(d.evidence.iter().any(|e| e == "ltr_of:int"));
        // without a homology-backed label for the partner, nothing changes
        let d = apply_neighbours(unknown(), adj.get("ltr"), &HashMap::new(), &p);
        assert_eq!(d.label, "Unknown");
    }

    #[test]
    fn distant_or_rare_neighbours_are_not_evidence() {
        let mut by_family: HashMap<String, Vec<RmHit>> = HashMap::new();
        for i in 0..10u64 {
            let c = format!("c{i}");
            by_family.entry("u".into()).or_default().push(hit(&c, 1000, 1500, false, "u"));
            // 200 bp away: not adjacent
            by_family.entry("dna".into()).or_default().push(hit(&c, 1700, 3000, false, "dna"));
        }
        // adjacent in 2 of 10 copies only
        for i in 0..2u64 {
            let c = format!("c{i}");
            by_family.entry("line".into()).or_default().push(hit(&c, 500, 990, false, "line"));
        }
        let p = NeighbourParams::default();
        let adj = adjacency(&by_family, &p);
        assert!(!adj["u"].partners.contains_key("dna"));
        let labels = HashMap::from([
            ("dna".to_owned(), "DNA/hAT".to_owned()),
            ("line".to_owned(), "LINE/L1".to_owned()),
        ]);
        let d = apply_neighbours(unknown(), adj.get("u"), &labels, &p);
        assert_eq!(d.label, "Unknown");
    }
}
