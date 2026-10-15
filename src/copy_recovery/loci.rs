use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use super::{Hit, Strand};

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct Locus {
    pub family_id: String,
    pub locus_id: String,
    pub subject_id: String,
    pub start: u64,
    pub end: u64,
    pub extended_start: u64,
    pub extended_end: u64,
    pub strand: Strand,
    pub n_fragments: usize,
    pub best_bitscore: f64,
    pub combined_query_coverage: f64,
    pub percent_identity: f64,
    pub representative_score: f64,
    pub selected: bool,
    pub exclusion_reason: Option<String>,
    pub orientation_method: String,
    pub warnings: Vec<String>,
    #[serde(skip)]
    pub fragments: Vec<Hit>,
}

/// Query overlap (bp) tolerated between consecutive fragments of one locus.
const MAX_QUERY_OVERLAP: u64 = 50;

fn coherent(previous: &Hit, next: &Hit, max_gap: u64) -> bool {
    if previous.family_id != next.family_id
        || previous.subject_id != next.subject_id
        || previous.strand != next.strand
        || next.subject_start.saturating_sub(previous.subject_end) > max_gap
    {
        return false;
    }
    // Fragments of one interrupted copy cover successive, essentially
    // non-overlapping query ranges. Two separate copies that both start near
    // the query 5' end must not be merged just because they are close.
    match previous.strand {
        Strand::Plus => next.query_start + MAX_QUERY_OVERLAP > previous.query_end,
        Strand::Minus => next.query_end < previous.query_start + MAX_QUERY_OVERLAP,
    }
}

pub fn merge_hits(hits: &[Hit], max_gap: u64) -> Vec<Locus> {
    let mut groups: BTreeMap<(String, String, Strand), Vec<Hit>> = BTreeMap::new();
    for hit in hits {
        groups
            .entry((hit.family_id.clone(), hit.subject_id.clone(), hit.strand))
            .or_default()
            .push(hit.clone());
    }
    let mut loci = Vec::new();
    for ((family, subject, strand), mut group) in groups {
        group.sort_by_key(|hit| (hit.subject_start, hit.subject_end, hit.query_start));
        let mut fragments: Vec<Hit> = Vec::new();
        for hit in group {
            if !fragments.is_empty()
                && !coherent(fragments.last().expect("nonempty"), &hit, max_gap)
            {
                loci.push(make_locus(
                    &family,
                    &subject,
                    strand,
                    std::mem::take(&mut fragments),
                ));
            }
            fragments.push(hit);
        }
        if !fragments.is_empty() {
            loci.push(make_locus(&family, &subject, strand, fragments));
        }
    }
    assign_ids(&mut loci);
    loci
}

fn make_locus(family: &str, subject: &str, strand: Strand, fragments: Vec<Hit>) -> Locus {
    let start = fragments.iter().map(|h| h.subject_start).min().unwrap() - 1;
    let end = fragments.iter().map(|h| h.subject_end).max().unwrap();
    let best_bitscore = fragments.iter().map(|h| h.bitscore).fold(0.0, f64::max);
    let query_length = fragments.first().map_or(1, |h| h.query_length.max(1));
    let mut query_intervals = fragments
        .iter()
        .map(|h| (h.query_start - 1, h.query_end))
        .collect::<Vec<_>>();
    query_intervals.sort_unstable();
    let mut covered = 0;
    let mut last_end = 0;
    for (start, end) in query_intervals {
        if end > last_end {
            covered += end - start.max(last_end);
            last_end = end;
        }
    }
    let coverage = (covered as f64 / query_length as f64).min(1.0);
    let identity =
        fragments.iter().map(|h| h.percent_identity).sum::<f64>() / fragments.len() as f64;
    let fragmentation_penalty = (fragments.len().saturating_sub(1) as f64) * 2.0;
    let score = coverage * 50.0 + identity * 0.3 + best_bitscore.log10().max(0.0) * 10.0
        - fragmentation_penalty;
    Locus {
        family_id: family.into(),
        locus_id: String::new(),
        subject_id: subject.into(),
        start,
        end,
        extended_start: start,
        extended_end: end,
        strand,
        n_fragments: fragments.len(),
        best_bitscore,
        combined_query_coverage: coverage,
        percent_identity: identity,
        representative_score: score,
        selected: false,
        exclusion_reason: None,
        orientation_method: "blast-strand".into(),
        warnings: Vec::new(),
        fragments,
    }
}

fn assign_ids(loci: &mut [Locus]) {
    loci.sort_by(|a, b| {
        (&a.family_id, &a.subject_id, a.start, a.end, a.strand).cmp(&(
            &b.family_id,
            &b.subject_id,
            b.start,
            b.end,
            b.strand,
        ))
    });
    for (index, locus) in loci.iter_mut().enumerate() {
        locus.locus_id = format!("locus-{:06}", index + 1);
    }
}

/// Collapse overlapping loci of the same family on the same subject.
///
/// Every input locus is retained: the best-scoring locus in an overlapping
/// cluster stays eligible and the others are marked `duplicate_of:<locus_id>`
/// so rejected evidence remains auditable. Loci that were already excluded
/// (for example as ambiguous-strand) never act as a cluster representative.
/// Locus IDs assigned by [`merge_hits`] are preserved.
pub fn collapse_overlapping_loci(mut loci: Vec<Locus>) -> Vec<Locus> {
    loci.sort_by(|a, b| {
        (&a.family_id, &a.subject_id, a.start, &a.locus_id).cmp(&(
            &b.family_id,
            &b.subject_id,
            b.start,
            &b.locus_id,
        ))
    });
    let mut representative: Option<usize> = None;
    for index in 0..loci.len() {
        if loci[index].exclusion_reason.is_some() {
            continue;
        }
        let Some(rep) = representative else {
            representative = Some(index);
            continue;
        };
        let (previous, locus) = (&loci[rep], &loci[index]);
        let overlap = previous.family_id == locus.family_id
            && previous.subject_id == locus.subject_id
            && locus.start < previous.end
            && (previous.end.min(locus.end) - locus.start) * 2
                >= (previous.end - previous.start).min(locus.end - locus.start);
        if !overlap {
            representative = Some(index);
            continue;
        }
        let (loser, winner) = if locus.representative_score > previous.representative_score {
            (rep, index)
        } else {
            (index, rep)
        };
        let reason = format!("duplicate_of:{}", loci[winner].locus_id);
        loci[loser].selected = false;
        loci[loser].exclusion_reason = Some(reason);
        representative = Some(winner);
    }
    loci.sort_by(|a, b| a.locus_id.cmp(&b.locus_id));
    loci
}

/// Exclude overlapping opposite-strand loci whose coverage is too weak to
/// decide the copy orientation. When one covers at least twice as much of the
/// query, it decides: the other is the shadow of a terminal inverted repeat (or
/// another palindromic stretch) matching the query backwards, and only it is
/// excluded, as `opposite_strand_shadow:<locus_id>`.
pub fn mark_ambiguous_strand_loci(loci: &mut [Locus], threshold: f64) {
    for i in 0..loci.len() {
        for j in (i + 1)..loci.len() {
            let (ci, cj) = (
                loci[i].combined_query_coverage,
                loci[j].combined_query_coverage,
            );
            if loci[i].family_id == loci[j].family_id
                && loci[i].subject_id == loci[j].subject_id
                && loci[i].strand != loci[j].strand
                && loci[i].start < loci[j].end
                && loci[j].start < loci[i].end
                && ci.min(cj) < threshold
            {
                let shadow = if ci >= 2.0 * cj {
                    Some((j, i))
                } else if cj >= 2.0 * ci {
                    Some((i, j))
                } else {
                    None
                };
                if let Some((weak, strong)) = shadow {
                    let reason = format!("opposite_strand_shadow:{}", loci[strong].locus_id);
                    loci[weak].selected = false;
                    loci[weak].exclusion_reason = Some(reason);
                    continue;
                }
                for k in [i, j] {
                    loci[k].selected = false;
                    loci[k].exclusion_reason = Some("ambiguous_strand".into());
                }
            }
        }
    }
}

/// Select up to `limit` eligible loci by descending score with deterministic
/// tie-breaking. Loci that already carry an exclusion reason are never selected.
pub fn rank_loci(loci: &mut [Locus], limit: usize, diversity_filter: bool) {
    let mut order = (0..loci.len())
        .filter(|&index| loci[index].exclusion_reason.is_none())
        .collect::<Vec<_>>();
    order.sort_by(|&a, &b| {
        loci[b]
            .representative_score
            .total_cmp(&loci[a].representative_score)
            .then_with(|| loci[a].subject_id.cmp(&loci[b].subject_id))
            .then_with(|| loci[a].start.cmp(&loci[b].start))
            .then_with(|| loci[a].locus_id.cmp(&loci[b].locus_id))
    });
    let mut selected: Vec<usize> = Vec::new();
    for index in order {
        if selected.len() == limit {
            loci[index].exclusion_reason = Some("selection_limit".into());
            continue;
        }
        if diversity_filter
            && selected.iter().any(|&other| {
                loci[other].subject_id == loci[index].subject_id
                    && loci[other].start.abs_diff(loci[index].start) < 100
            })
        {
            loci[index].exclusion_reason = Some("redundant_locus".into());
        } else {
            loci[index].selected = true;
            selected.push(index);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn hit(qs: u64, qe: u64, ss: u64, se: u64, strand: Strand) -> Hit {
        Hit {
            family_id: "f".into(),
            query_id: "f".into(),
            subject_id: "chr".into(),
            query_start: qs,
            query_end: qe,
            query_length: 100,
            subject_start: ss.min(se),
            subject_end: ss.max(se),
            strand,
            alignment_length: qe - qs + 1,
            matches: qe - qs + 1,
            mismatches: 0,
            gap_opens: 0,
            percent_identity: 100.0,
            query_coverage: (qe - qs + 1) as f64 / 100.0,
            bitscore: 100.0,
            evalue: 1e-20,
        }
    }
    #[test]
    fn merges_coherent_fragments() {
        assert_eq!(
            merge_hits(
                &[
                    hit(1, 40, 10, 49, Strand::Plus),
                    hit(41, 80, 55, 94, Strand::Plus)
                ],
                10
            )
            .len(),
            1
        );
    }
    #[test]
    fn rejects_incoherent_query_order() {
        assert_eq!(
            merge_hits(
                &[
                    hit(60, 90, 10, 40, Strand::Plus),
                    hit(1, 30, 45, 74, Strand::Plus)
                ],
                10
            )
            .len(),
            2
        );
    }
    #[test]
    fn minus_order_is_coherent_when_query_decreases() {
        assert_eq!(
            merge_hits(
                &[
                    hit(61, 100, 10, 49, Strand::Minus),
                    hit(21, 60, 55, 94, Strand::Minus)
                ],
                10
            )
            .len(),
            1
        );
    }
    #[test]
    fn collapse_retains_duplicates_with_reason() {
        let mut low = hit(1, 50, 10, 59, Strand::Plus);
        low.bitscore = 10.0;
        let loci = merge_hits(&[low], 10);
        let mut high = merge_hits(&[hit(1, 100, 20, 119, Strand::Plus)], 10);
        high[0].locus_id = "locus-000002".into();
        let all = collapse_overlapping_loci(loci.into_iter().chain(high).collect());
        assert_eq!(all.len(), 2);
        assert_eq!(
            all[0].exclusion_reason.as_deref(),
            Some("duplicate_of:locus-000002")
        );
        assert!(all[1].exclusion_reason.is_none());
    }
    #[test]
    fn ranking_never_selects_excluded_loci() {
        let mut loci = merge_hits(
            &[
                hit(1, 30, 10, 39, Strand::Plus),
                hit(1, 30, 20, 49, Strand::Minus),
            ],
            10,
        );
        mark_ambiguous_strand_loci(&mut loci, 0.8);
        rank_loci(&mut loci, 10, false);
        assert!(loci.iter().all(|l| !l.selected));
        assert!(loci
            .iter()
            .all(|l| l.exclusion_reason.as_deref() == Some("ambiguous_strand")));
    }
    #[test]
    fn a_full_length_copy_outweighs_its_inverted_repeat_shadow() {
        // a TIR element: the whole query matches the copy on the plus strand,
        // and the 5' TIR also matches its 3' end backwards
        let mut loci = merge_hits(
            &[
                hit(1, 100, 1000, 1099, Strand::Plus),
                hit(1, 20, 1080, 1099, Strand::Minus),
            ],
            10,
        );
        mark_ambiguous_strand_loci(&mut loci, 0.8);
        rank_loci(&mut loci, 10, false);
        let plus = loci.iter().find(|l| l.strand == Strand::Plus).unwrap();
        let minus = loci.iter().find(|l| l.strand == Strand::Minus).unwrap();
        assert!(plus.selected && plus.exclusion_reason.is_none());
        assert_eq!(
            minus.exclusion_reason.as_deref(),
            Some(format!("opposite_strand_shadow:{}", plus.locus_id).as_str())
        );
    }
    #[test]
    fn ranking_limit_and_ties_are_deterministic() {
        let mut loci = merge_hits(
            &[
                hit(1, 100, 1000, 1099, Strand::Plus),
                hit(1, 100, 10, 109, Strand::Plus),
                hit(1, 100, 5000, 5099, Strand::Plus),
            ],
            10,
        );
        rank_loci(&mut loci, 2, false);
        let chosen = loci
            .iter()
            .filter(|l| l.selected)
            .map(|l| l.start)
            .collect::<Vec<_>>();
        assert_eq!(chosen, vec![9, 999]);
        assert_eq!(loci[2].exclusion_reason.as_deref(), Some("selection_limit"));
    }
    #[test]
    fn does_not_merge_neighbouring_full_and_partial_copies() {
        // full-length copy then a 5' truncated copy 8 kb downstream
        assert_eq!(
            merge_hits(
                &[
                    hit(1, 100, 1000, 1099, Strand::Plus),
                    hit(1, 60, 9000, 9059, Strand::Plus)
                ],
                10_000
            )
            .len(),
            2
        );
        assert_eq!(
            merge_hits(
                &[
                    hit(1, 100, 1000, 1099, Strand::Minus),
                    hit(40, 100, 9000, 9060, Strand::Minus)
                ],
                10_000
            )
            .len(),
            2
        );
    }
}
