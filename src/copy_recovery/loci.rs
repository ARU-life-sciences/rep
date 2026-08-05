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

fn coherent(previous: &Hit, next: &Hit, max_gap: u64) -> bool {
    if previous.family_id != next.family_id
        || previous.subject_id != next.subject_id
        || previous.strand != next.strand
        || next.subject_start.saturating_sub(previous.subject_end) > max_gap
    {
        return false;
    }
    match previous.strand {
        Strand::Plus => next.query_start >= previous.query_start,
        Strand::Minus => next.query_end <= previous.query_end,
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

pub fn collapse_overlapping_loci(mut loci: Vec<Locus>) -> Vec<Locus> {
    loci.sort_by(|a, b| {
        (&a.family_id, &a.subject_id, a.start).cmp(&(&b.family_id, &b.subject_id, b.start))
    });
    let mut result: Vec<Locus> = Vec::new();
    for mut locus in loci {
        if let Some(previous) = result.last_mut() {
            let overlap = previous.family_id == locus.family_id
                && previous.subject_id == locus.subject_id
                && locus.start < previous.end
                && (previous.end.min(locus.end) - locus.start) * 2
                    >= (previous.end - previous.start).min(locus.end - locus.start);
            if overlap {
                if locus.representative_score > previous.representative_score {
                    previous.exclusion_reason = Some(format!("collapsed_into:{}", locus.locus_id));
                    std::mem::swap(previous, &mut locus);
                }
                locus.selected = false;
                locus.exclusion_reason = Some(format!("duplicate_of:{}", previous.locus_id));
                continue;
            }
        }
        result.push(locus);
    }
    assign_ids(&mut result);
    result
}

pub fn mark_ambiguous_strand_loci(loci: &mut [Locus], threshold: f64) {
    for i in 0..loci.len() {
        for j in (i + 1)..loci.len() {
            if loci[i].family_id == loci[j].family_id
                && loci[i].subject_id == loci[j].subject_id
                && loci[i].strand != loci[j].strand
                && loci[i].start < loci[j].end
                && loci[j].start < loci[i].end
                && loci[i]
                    .combined_query_coverage
                    .min(loci[j].combined_query_coverage)
                    < threshold
            {
                loci[i].selected = false;
                loci[j].selected = false;
                loci[i].exclusion_reason = Some("ambiguous_strand".into());
                loci[j].exclusion_reason = Some("ambiguous_strand".into());
            }
        }
    }
}

pub fn rank_loci(loci: &mut [Locus], limit: usize, diversity_filter: bool) {
    let mut order = (0..loci.len()).collect::<Vec<_>>();
    order.sort_by(|&a, &b| {
        loci[b]
            .representative_score
            .total_cmp(&loci[a].representative_score)
            .then_with(|| loci[a].subject_id.cmp(&loci[b].subject_id))
            .then_with(|| loci[a].start.cmp(&loci[b].start))
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
}
