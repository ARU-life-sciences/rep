//! Repeats found from the sequence itself: tandem runs and low complexity
//! in unannotated stretches, and terminal-repeat pairs (LTR- or TIR-like)
//! that may mark elements no reference knows.

use crate::classify::{composition, tandem};

/// A tandem or low-complexity run (0-based half-open, relative to the
/// scanned slice).
#[derive(Clone, Debug, PartialEq)]
pub struct Run {
    pub start: usize,
    pub end: usize,
    /// Repeat unit length; 0 for low complexity.
    pub period: usize,
    pub identity: f64,
}

/// Microsatellite periods scanned base by base.
const MICRO_MAX_PERIOD: usize = 10;
/// Alignment scores of a base against the base one period on, and the
/// smallest run score reported. More permissive than TRF's 2, 7 and 50:
/// RepeatMasker reports simple repeats of ~40 bp at 5–25% divergence, which
/// those settings miss. Random sequence still gives almost no calls (tested).
const MICRO_MATCH: i64 = 2;
const MICRO_MISMATCH: i64 = -5;
const MICRO_MIN_SCORE: i64 = 30;

/// Satellites with units of 11–200 bp, by windows holding at least three
/// units: (window, step, min identity). Edges are accurate to the step.
const SATELLITE_SCALE: (usize, usize, f64) = (600, 200, 0.75);

/// Tandem runs in `seq`: microsatellites (units of 1–10 bp) found base by
/// base, then satellites (units of 11–200 bp) found by windows.
pub fn tandem_runs(seq: &[u8]) -> Vec<Run> {
    let mut runs = microsatellites(seq);
    runs.extend(satellites(seq));
    // microsatellites go first, so overlaps resolve in their favour downstream
    runs
}

/// Microsatellites: for each period, maximal-scoring stretches where bases
/// match the base one period downstream (match +2, mismatch or N -5, at
/// least 30). Overlapping calls keep the best score, then the shortest period,
/// so a (CA)n run is not reported again as (CACA)n.
pub fn microsatellites(seq: &[u8]) -> Vec<Run> {
    let mut calls: Vec<(i64, Run)> = Vec::new();
    for p in 1..=MICRO_MAX_PERIOD.min(seq.len().saturating_sub(1)) {
        let (mut score, mut start) = (0i64, 0usize);
        let (mut best, mut best_end, mut matches) = (0i64, 0usize, 0usize);
        let mut best_matches = 0;
        let mut flush = |best: i64, start: usize, best_end: usize, m: usize| {
            if best >= MICRO_MIN_SCORE {
                let len = best_end - start;
                calls.push((
                    best,
                    Run {
                        start,
                        end: best_end + p,
                        period: p,
                        identity: m as f64 / len as f64,
                    },
                ));
            }
        };
        for i in 0..seq.len() - p {
            let same = seq[i] == seq[i + p] && seq[i] != b'N';
            score += if same { MICRO_MATCH } else { MICRO_MISMATCH };
            matches += usize::from(same);
            if score > best {
                (best, best_end, best_matches) = (score, i + 1, matches);
            }
            if score <= 0 {
                flush(best, start, best_end, best_matches);
                (score, start, best, matches) = (0, i + 1, 0, 0);
            }
        }
        flush(best, start, best_end, best_matches);
    }
    calls.sort_by(|a, b| {
        b.0.cmp(&a.0)
            .then(a.1.period.cmp(&b.1.period))
            .then(a.1.start.cmp(&b.1.start))
    });
    let mut kept: Vec<Run> = Vec::new();
    for (_, r) in calls {
        if kept.iter().all(|k| r.end <= k.start || r.start >= k.end) {
            kept.push(r);
        }
    }
    kept.sort_by_key(|r| r.start);
    kept
}

/// Satellite runs, merging overlapping windows that call the same period
/// (within 10%).
fn satellites(seq: &[u8]) -> Vec<Run> {
    let (window, step, min_identity) = SATELLITE_SCALE;
    let mut runs = Vec::new();
    let mut current: Option<(Run, usize)> = None;
    let mut start = 0;
    while start + window <= seq.len() {
        let w = &seq[start..start + window];
        let call = (composition(w).n_fraction <= 0.1)
            .then(|| tandem(w, 200, window, min_identity))
            .flatten()
            .filter(|t| t.period > MICRO_MAX_PERIOD);
        match (call, current.as_mut()) {
            (Some(t), Some((run, windows)))
                if start <= run.end && t.period.abs_diff(run.period) * 10 <= run.period =>
            {
                run.end = start + window;
                run.identity += t.identity;
                *windows += 1;
            }
            (call, _) => {
                if let Some((run, windows)) = current.take() {
                    runs.push(Run {
                        identity: run.identity / windows as f64,
                        ..run
                    });
                }
                current = call.map(|t| {
                    (
                        Run {
                            start,
                            end: start + window,
                            period: t.period,
                            identity: t.identity,
                        },
                        1,
                    )
                });
            }
        }
        start += step;
    }
    if let Some((run, windows)) = current {
        runs.push(Run {
            identity: run.identity / windows as f64,
            ..run
        });
    }
    runs
}

/// Low-complexity runs (base-composition bias without periodicity).
pub fn low_complexity_runs(seq: &[u8]) -> Vec<Run> {
    let (window, step) = (40, 10);
    let mut runs: Vec<Run> = Vec::new();
    let mut start = 0;
    while start + window <= seq.len() {
        let c = composition(&seq[start..start + window]);
        if c.low_complexity && c.n_fraction <= 0.1 {
            match runs.last_mut() {
                Some(r) if start <= r.end => r.end = start + window,
                _ => runs.push(Run {
                    start,
                    end: start + window,
                    period: 0,
                    identity: 0.0,
                }),
            }
        }
        start += step;
    }
    runs
}

/// One self-alignment of an unannotated stretch (1-based inclusive, BLASTN
/// order: `s_start > s_end` on the minus strand).
#[derive(Clone, Debug)]
pub struct SelfHit {
    pub q_start: u64,
    pub q_end: u64,
    pub s_start: u64,
    pub s_end: u64,
    pub identity: f64,
}

/// A candidate element bounded by a terminal-repeat pair, in stretch
/// coordinates (0-based half-open).
#[derive(Clone, Debug, PartialEq)]
pub struct Candidate {
    pub kind: &'static str,
    pub start: u64,
    pub end: u64,
    pub repeat_length: u64,
    pub identity: f64,
}

pub const LTR_PAIR: &str = "candidate_LTR_pair";
pub const TIR_PAIR: &str = "candidate_TIR_pair";

/// Terminal-repeat pairs among a stretch's self-alignments:
/// - direct repeats of 100–3,000 bp, at least 80% identical, 1–25 kb apart
///   with at least 500 bp between them (LTR retrotransposon-like);
/// - inverted repeats of 20–1,000 bp, at least 80% identical, bounding
///   80 bp–20 kb with at least 50 bp between them (TIR transposon-like).
///
/// Overlapping candidates keep the longest, most identical pair.
pub fn terminal_pairs(hits: &[SelfHit]) -> Vec<Candidate> {
    let mut all = Vec::new();
    for h in hits {
        let len = h.q_end + 1 - h.q_start;
        if h.identity < 80.0 {
            continue;
        }
        if h.s_start <= h.s_end {
            // direct: report each pair once, with the query copy first
            if h.s_start <= h.q_end {
                continue;
            }
            let inner = h.s_start - h.q_end - 1;
            let span = h.s_end - h.q_start + 1;
            if (100..=3000).contains(&len) && inner >= 500 && (1000..=25_000).contains(&span) {
                all.push(Candidate {
                    kind: LTR_PAIR,
                    start: h.q_start - 1,
                    end: h.s_end,
                    repeat_length: len,
                    identity: h.identity,
                });
            }
        } else {
            // inverted: the subject copy (s_end..s_start) lies after the query copy
            if h.s_end <= h.q_end {
                continue;
            }
            let inner = h.s_end - h.q_end - 1;
            let span = h.s_start - h.q_start + 1;
            if (20..=1000).contains(&len) && inner >= 50 && (80..=20_000).contains(&span) {
                all.push(Candidate {
                    kind: TIR_PAIR,
                    start: h.q_start - 1,
                    end: h.s_start,
                    repeat_length: len,
                    identity: h.identity,
                });
            }
        }
    }
    all.sort_by(|a, b| {
        (b.repeat_length as f64 * b.identity)
            .total_cmp(&(a.repeat_length as f64 * a.identity))
            .then(a.start.cmp(&b.start))
    });
    let mut kept: Vec<Candidate> = Vec::new();
    for c in all {
        if kept.iter().all(|k| c.end <= k.start || c.start >= k.end) {
            kept.push(c);
        }
    }
    kept.sort_by_key(|c| c.start);
    kept
}

/// Length of the longest target-site duplication (`min..=max` bp) flanking
/// `[start, end)` in `seq`, if any.
pub fn tsd(seq: &[u8], start: usize, end: usize, min: usize, max: usize) -> Option<usize> {
    (min..=max).rev().find(|&k| {
        start >= k
            && end + k <= seq.len()
            && seq[start - k..start] == seq[end..end + k]
            && !seq[end..end + k].contains(&b'N')
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 62) as usize]
            })
            .collect()
    }

    #[test]
    fn microsatellites_and_satellites_are_found_in_random_sequence() {
        let mut seq = random(2000, 1);
        seq.extend(b"CA".repeat(60));
        seq.extend(random(2000, 2));
        let unit = random(150, 3);
        for _ in 0..8 {
            seq.extend(&unit);
        }
        seq.extend(random(2000, 4));
        let runs = tandem_runs(&seq);
        let micro = runs.iter().find(|r| r.period == 2).expect("CA run");
        assert!(micro.start >= 1960 && micro.end <= 2140, "{micro:?}");
        let sat = runs
            .iter()
            .find(|r| r.period == 150)
            .expect("150 bp satellite");
        assert!(sat.start >= 3900 && sat.end <= 5540, "{sat:?}");
        assert_eq!(runs.len(), 2, "{runs:?}");
        // a 45 bp (GAA)n with one substitution, as TRF would report it
        let mut short = random(300, 7);
        short.extend(b"GAAGAAGAAGAAGTAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAA");
        short.extend(random(300, 8));
        let m = microsatellites(&short);
        assert_eq!(m.len(), 1, "{m:?}");
        assert_eq!(m[0].period, 3);
        // random bases next to the run may happen to continue it
        assert!((297..=300).contains(&m[0].start), "{m:?}");
        assert!(m[0].end >= 343 && m[0].end <= 347, "{m:?}");
        assert!(low_complexity_runs(&random(3000, 5)).is_empty());
        let mut biased = random(500, 6);
        biased.extend(b"AAAGAAAAGAAAAAGAAAGAAAAAGAAAAGAAAAAAGAAAAGAAAAAAGAAAAAGAAAAGAAAA");
        assert_eq!(low_complexity_runs(&biased).len(), 1);
    }

    #[test]
    fn random_sequence_has_almost_no_simple_repeats() {
        let seq = random(1_000_000, 11);
        let bp = |runs: Vec<Run>| runs.iter().map(|r| r.end - r.start).sum::<usize>();
        let (micro, low) = (bp(microsatellites(&seq)), bp(low_complexity_runs(&seq)));
        assert!(micro < 500, "{micro} bp of microsatellites in 1 Mb");
        assert_eq!(low, 0);
    }

    #[test]
    fn terminal_pairs_are_ltr_or_tir_like() {
        let hit = |qs, qe, ss, se, id| SelfHit {
            q_start: qs,
            q_end: qe,
            s_start: ss,
            s_end: se,
            identity: id,
        };
        let c = terminal_pairs(&[
            // 400 bp LTRs, 5 kb element
            hit(1001, 1400, 5601, 6000, 95.0),
            // the same pair reported the other way round
            hit(5601, 6000, 1001, 1400, 95.0),
            // 30 bp TIRs bounding 1 kb
            hit(8001, 8030, 9000, 8971, 90.0),
            // a tandem duplication: no inner region
            hit(12001, 12400, 12401, 12800, 99.0),
            // too diverged
            hit(15001, 15100, 15800, 15701, 70.0),
        ]);
        assert_eq!(
            c,
            vec![
                Candidate {
                    kind: LTR_PAIR,
                    start: 1000,
                    end: 6000,
                    repeat_length: 400,
                    identity: 95.0
                },
                Candidate {
                    kind: TIR_PAIR,
                    start: 8000,
                    end: 9000,
                    repeat_length: 30,
                    identity: 90.0
                },
            ]
        );
        let s = b"GGGGCATGAAAAAAAAAACATGTTTT";
        assert_eq!(tsd(s, 8, 18, 2, 10), Some(4));
        assert_eq!(tsd(s, 9, 18, 2, 10), None);
    }
}
