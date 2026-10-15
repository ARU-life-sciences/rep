//! Overlap resolution and fragment chaining.
//!
//! Coordinates are 0-based half-open on the forward strand of the input
//! sequence; reference coordinates are 1-based inclusive, as BLAST reports
//! them (amino acids for protein references).

use std::collections::{BTreeMap, HashMap};

use super::{Evidence, Fragment};

/// Bases already annotated on one sequence, as disjoint intervals.
#[derive(Clone, Debug, Default)]
pub struct Occupancy(BTreeMap<u64, u64>);

impl Occupancy {
    /// The parts of `[start, end)` not yet annotated.
    pub fn free(&self, start: u64, end: u64) -> Vec<(u64, u64)> {
        let mut pieces = Vec::new();
        let mut at = start;
        // the interval starting before `start` may still cover it
        let first = self.0.range(..start).next_back().map(|(&s, &e)| (s, e));
        for (s, e) in first
            .into_iter()
            .chain(self.0.range(start..end).map(|(&s, &e)| (s, e)))
        {
            if e <= at {
                continue;
            }
            if s > at {
                pieces.push((at, s.min(end)));
            }
            at = at.max(e);
            if at >= end {
                break;
            }
        }
        if at < end {
            pieces.push((at, end));
        }
        pieces
    }

    /// Mark `[start, end)`, which must be free, as annotated.
    pub fn insert(&mut self, start: u64, end: u64) {
        self.0.insert(start, end);
    }

    /// Unannotated stretches of a sequence of length `len`, at least `min` long.
    pub fn gaps(&self, len: u64, min: u64) -> Vec<(u64, u64)> {
        self.free(0, len)
            .into_iter()
            .filter(|(s, e)| e - s >= min)
            .collect()
    }
}

/// Give each base to the best-scoring hit covering it, in the order given
/// (callers sort by score). A hit keeps its unclaimed pieces of at least
/// `min_length`; reference coordinates of a piece are interpolated.
pub fn claim(
    hits: Vec<Fragment>,
    occupied: &mut HashMap<usize, Occupancy>,
    min_length: u64,
) -> Vec<Fragment> {
    let mut kept = Vec::new();
    for hit in hits {
        let occ = occupied.entry(hit.seq).or_default();
        for (s, e) in occ.free(hit.start, hit.end) {
            if e - s < min_length {
                continue;
            }
            occ.insert(s, e);
            kept.push(trim(&hit, s, e));
        }
    }
    kept
}

/// The part `[s, e)` of a hit, with its reference stretch scaled to match.
fn trim(hit: &Fragment, s: u64, e: u64) -> Fragment {
    let mut piece = hit.clone();
    piece.start = s;
    piece.end = e;
    if hit.ref_length == 0 || (s == hit.start && e == hit.end) {
        return piece;
    }
    let span = (hit.end - hit.start) as f64;
    let ref_span = (hit.ref_end + 1 - hit.ref_start) as f64;
    let at = |q: u64| ((q - hit.start) as f64 / span * ref_span).round() as u64;
    let (a, b) = (at(s), at(e));
    if hit.strand == '+' {
        piece.ref_start = hit.ref_start + a;
        piece.ref_end = (hit.ref_start + b).saturating_sub(1).max(piece.ref_start);
    } else {
        piece.ref_end = hit.ref_end.saturating_sub(a);
        piece.ref_start = (hit.ref_end + 1).saturating_sub(b).min(piece.ref_end);
    }
    piece
}

/// Number fragments into elements: fragments of the same reference on the
/// same strand join one element when they follow each other along the
/// reference and are at most `max_gap` apart once other annotated fragments
/// between them (insertions) are discounted. Fragments must be sorted by
/// sequence and start and must not overlap. Returns the number of elements.
pub fn chain(fragments: &mut [Fragment], max_gap: u64) -> usize {
    // annotated bases before each fragment, per sequence
    let mut before = Vec::with_capacity(fragments.len());
    let mut total = 0;
    let mut seq = usize::MAX;
    for f in fragments.iter() {
        if f.seq != seq {
            seq = f.seq;
            total = 0;
        }
        before.push(total);
        total += f.end - f.start;
    }
    let mut open: HashMap<(usize, String, char), usize> = HashMap::new();
    let mut elements = 0;
    for i in 0..fragments.len() {
        let f = &fragments[i];
        let homology = matches!(f.evidence, Evidence::Nucleotide | Evidence::Protein);
        let key = (f.seq, f.reference.clone(), f.strand);
        let joined = homology
            .then(|| open.get(&key).copied())
            .flatten()
            .filter(|&j| {
                let last = &fragments[j];
                let between = before[i] - before[j] - (last.end - last.start);
                let gap = f.start.saturating_sub(last.end).saturating_sub(between);
                let in_order = if f.strand == '+' {
                    f.ref_start >= last.ref_start
                } else {
                    f.ref_end <= last.ref_end
                };
                gap <= max_gap && in_order
            });
        let element = match joined {
            Some(j) => fragments[j].element,
            None => {
                elements += 1;
                elements
            }
        };
        fragments[i].element = element;
        if homology {
            open.insert(key, i);
        }
    }
    elements
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frag(start: u64, end: u64, reference: &str, strand: char, rs: u64, re: u64) -> Fragment {
        Fragment {
            seq: 0,
            start,
            end,
            strand,
            reference: reference.into(),
            ref_start: rs,
            ref_end: re,
            ref_length: 5000,
            class: "LTR/Gypsy".into(),
            evidence: Evidence::Nucleotide,
            ..Default::default()
        }
    }

    #[test]
    fn free_pieces_skip_annotated_bases() {
        let mut o = Occupancy::default();
        o.insert(10, 20);
        o.insert(30, 40);
        assert_eq!(o.free(0, 50), vec![(0, 10), (20, 30), (40, 50)]);
        assert_eq!(o.free(15, 35), vec![(20, 30)]);
        assert_eq!(o.free(12, 18), vec![]);
        assert_eq!(o.gaps(50, 11), vec![]);
        assert_eq!(o.gaps(45, 10), vec![(0, 10), (20, 30)]);
    }

    #[test]
    fn best_hits_claim_bases_and_losers_keep_their_rest() {
        let mut occ = HashMap::new();
        let best = frag(100, 200, "A", '+', 1, 100);
        let worse = frag(150, 350, "B", '-', 1, 200);
        let kept = claim(vec![best, worse], &mut occ, 30);
        assert_eq!(kept.len(), 2);
        let b = &kept[1];
        assert_eq!((b.start, b.end), (200, 350));
        // minus strand: the query end maps to the reference start
        assert_eq!((b.ref_start, b.ref_end), (1, 150));
        // a piece below the minimum length is dropped
        let kept = claim(vec![frag(190, 215, "C", '+', 1, 25)], &mut occ, 30);
        assert!(kept.is_empty());
    }

    #[test]
    fn fragments_chain_across_insertions_in_reference_order() {
        let mut f = vec![
            frag(0, 1000, "Gypsy-1", '+', 1, 1000),
            // a 3 kb insertion of another element
            frag(1050, 4050, "L1-2", '-', 1, 3000),
            frag(4100, 5000, "Gypsy-1", '+', 1100, 2000),
            // the same reference again, but back at its start: a new copy
            frag(5100, 5500, "Gypsy-1", '+', 1, 400),
        ];
        assert_eq!(chain(&mut f, 200), 3);
        assert_eq!(
            f.iter().map(|x| x.element).collect::<Vec<_>>(),
            vec![1, 2, 1, 3]
        );
    }
}
