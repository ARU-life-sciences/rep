//! Evidence from recovered genomic copies (`rep recover-copies` output).
//!
//! For each family: copy number and completeness from `loci.tsv`, divergence
//! from the consensus, target-site duplications (TSDs) at the ends of
//! full-length copies (pinned by aligning the consensus termini), and whether the flanks of different copies resemble
//! each other (expected for segmental duplications or truncated consensi, not
//! for transposable elements).

use std::{collections::HashMap, path::Path};

use bio::{alignment::pairwise::Aligner, io::fasta};
use serde::Serialize;

use crate::Result;

#[derive(Clone, Debug, Serialize)]
pub struct CopyParams {
    /// Combined consensus coverage for a copy to count as full-length.
    pub full_length_coverage: f64,
    pub tsd_min: usize,
    pub tsd_max: usize,
    /// Shift (bp) into the flank used for the null TSD expectation.
    pub tsd_null_shift: usize,
    /// Minimum full-length copies needed to call a TSD.
    pub tsd_min_copies: usize,
    /// Flank length compared between copies.
    pub flank_window: usize,
}

impl Default for CopyParams {
    fn default() -> Self {
        Self {
            full_length_coverage: 0.8,
            tsd_min: 2,
            tsd_max: 15,
            tsd_null_shift: 40,
            tsd_min_copies: 3,
            flank_window: 50,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct Tsd {
    pub length: usize,
    /// Fraction of evaluated copies with this TSD length (±1 bp).
    pub support: f64,
    /// The same fraction at positions shifted into the flanks.
    pub null_support: f64,
    /// Most common TSD sequence at the modal length.
    pub motif: String,
    pub copies_evaluated: usize,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct CopyEvidence {
    /// Independent loci (duplicates and ambiguous-strand loci excluded).
    pub loci: usize,
    pub full_length: usize,
    /// Mean percent identity of loci to the consensus.
    pub mean_identity: f64,
    pub tsd: Option<Tsd>,
    /// Mean pairwise identity of copy flanks (0–1); `None` with < 2 copies.
    pub flank_similarity: Option<f64>,
    /// Genomic base pairs annotated as this family (RepeatMasker input only).
    pub genome_bp: u64,
    /// Helitron end and insertion-site signature on the copies.
    pub helitron_ends: Option<HelitronEnds>,
}

/// Read copy evidence for one family from its copy-recovery directory.
/// Returns `None` when the directory does not exist.
pub fn read_copy_evidence(
    dir: &Path,
    consensus: &[u8],
    params: &CopyParams,
) -> Result<Option<CopyEvidence>> {
    if !dir.join("loci.tsv").is_file() {
        return Ok(None);
    }
    let mut evidence = CopyEvidence::default();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(dir.join("loci.tsv"))?;
    let headers = reader.headers()?.clone();
    let col = |name: &str| headers.iter().position(|h| h == name);
    let (Some(cov_col), Some(id_col), Some(excl_col)) = (
        col("combined_query_coverage"),
        col("percent_identity"),
        col("exclusion_reason"),
    ) else {
        return Ok(None);
    };
    let mut identity_sum = 0.0;
    for row in reader.records() {
        let row = row?;
        let reason = row.get(excl_col).unwrap_or("");
        if reason.starts_with("duplicate_of")
            || reason.starts_with("opposite_strand_shadow")
            || reason == "ambiguous_strand"
        {
            continue;
        }
        evidence.loci += 1;
        identity_sum += row[id_col].parse::<f64>().unwrap_or(0.0);
        if row[cov_col].parse::<f64>().unwrap_or(0.0) >= params.full_length_coverage {
            evidence.full_length += 1;
        }
    }
    if evidence.loci > 0 {
        evidence.mean_identity = identity_sum / evidence.loci as f64;
    }

    let copies = oriented_copies(dir)?;
    let full = full_length_loci(dir, params)?;
    let complete = copies
        .iter()
        .filter(|c| full.contains_key(&c.locus_id))
        .collect::<Vec<_>>();
    evidence.tsd = call_tsd(&complete, consensus, params);
    evidence.helitron_ends = helitron_ends(&complete);
    evidence.flank_similarity = flank_similarity(&copies, params.flank_window);
    Ok(Some(evidence))
}

/// A selected copy in consensus orientation with its flanks.
#[derive(Clone, Debug, PartialEq)]
pub struct OrientedCopy {
    pub locus_id: String,
    pub seq: Vec<u8>,
    /// Length of the 5' (left, in consensus orientation) flank.
    pub left: usize,
    pub right: usize,
}

impl OrientedCopy {
    fn core_end(&self) -> usize {
        self.seq.len() - self.right
    }
}

fn oriented_copies(dir: &Path) -> Result<Vec<OrientedCopy>> {
    let fa = dir.join("copies.extended.fa");
    let tsv = dir.join("copies.tsv");
    if !fa.is_file() || !tsv.is_file() {
        return Ok(Vec::new());
    }
    let mut flanks: HashMap<String, (usize, usize)> = HashMap::new();
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(tsv)?;
    let headers = reader.headers()?.clone();
    let col = |name: &str| headers.iter().position(|h| h == name);
    let (Some(id), Some(left), Some(right), Some(rc)) = (
        col("locus_id"),
        col("left_flank_actual"),
        col("right_flank_actual"),
        col("reverse_complemented"),
    ) else {
        return Ok(Vec::new());
    };
    for row in reader.records() {
        let row = row?;
        let (l, r): (usize, usize) = (row[left].parse()?, row[right].parse()?);
        // flanks are reported in genomic orientation; minus-strand copies
        // were reverse-complemented, which swaps them
        let oriented = if &row[rc] == "true" { (r, l) } else { (l, r) };
        flanks.insert(row[id].to_owned(), oriented);
    }
    let mut copies = Vec::new();
    for record in fasta::Reader::from_file(fa)?.records() {
        let record = record?;
        let locus = record.id().trim_start_matches("locus_id=").to_owned();
        if let Some(&(left, right)) = flanks.get(&locus) {
            if left + right <= record.seq().len() {
                copies.push(OrientedCopy {
                    locus_id: locus,
                    seq: record.seq().to_ascii_uppercase(),
                    left,
                    right,
                });
            }
        }
    }
    Ok(copies)
}

fn full_length_loci(dir: &Path, params: &CopyParams) -> Result<HashMap<String, f64>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(dir.join("loci.tsv"))?;
    let headers = reader.headers()?.clone();
    let (Some(id), Some(cov)) = (
        headers.iter().position(|h| h == "locus_id"),
        headers.iter().position(|h| h == "combined_query_coverage"),
    ) else {
        return Ok(HashMap::new());
    };
    let mut full = HashMap::new();
    for row in reader.records() {
        let row = row?;
        let c: f64 = row[cov].parse().unwrap_or(0.0);
        if c >= params.full_length_coverage {
            full.insert(row[id].to_owned(), c);
        }
    }
    Ok(full)
}

/// Consensus bases aligned at each end of a copy to pin it.
const PIN_END: usize = 40;
/// Bases searched on either side of a reported copy boundary.
const PIN_WINDOW: usize = 60;
/// Consensus bases a pinning alignment may leave unaligned at the outer end.
const PIN_SLACK: usize = 5;

/// Element start and end in `c.seq`, from local alignments of the consensus
/// termini near the reported boundaries. Those come from BLAST or RepeatMasker
/// hits, which often stop a few bases inside divergent ends; a TSD has to be
/// sought exactly outside the element, where a few bases of slack would let
/// chance short matches in.
pub fn pin_ends(c: &OrientedCopy, consensus: &[u8]) -> Option<(usize, usize)> {
    let core = c.core_end().checked_sub(c.left)?;
    let k = PIN_END.min(consensus.len() / 2).min(core / 2);
    if k < 20 {
        return None;
    }
    let score = |a: u8, b: u8| {
        if a.eq_ignore_ascii_case(&b) {
            2i32
        } else {
            -3i32
        }
    };
    let mut aligner = Aligner::new(-5, -2, score);
    let window = |at: usize| {
        (
            at.saturating_sub(PIN_WINDOW),
            (at + PIN_WINDOW).min(c.seq.len()),
        )
    };
    // half the score of a perfect match
    let min_score = k as i32;
    let (a, b) = window(c.left);
    let five = aligner.local(&consensus[..k], &c.seq[a..b]);
    if five.score < min_score || five.xstart > PIN_SLACK {
        return None;
    }
    let start = (a + five.ystart).checked_sub(five.xstart)?;
    let (a, b) = window(c.core_end());
    let three = aligner.local(&consensus[consensus.len() - k..], &c.seq[a..b]);
    if three.score < min_score || k - three.xend > PIN_SLACK {
        return None;
    }
    let end = a + three.yend + (k - three.xend);
    (start < end && end <= c.seq.len()).then_some((start, end))
}

/// Longest duplication of the bases just before `start` just after `end`
/// (one mismatch allowed from 6 bp, as old TSDs have point mutations).
fn tsd_at(seq: &[u8], start: usize, end: usize, p: &CopyParams) -> Option<(usize, Vec<u8>)> {
    (p.tsd_min..=p.tsd_max).rev().find_map(|k| {
        if start < k || end + k > seq.len() {
            return None;
        }
        let left = &seq[start - k..start];
        let right = &seq[end..end + k];
        let mismatches = left.iter().zip(right).filter(|(a, b)| a != b).count();
        (mismatches <= usize::from(k >= 6) && !left.contains(&b'N') && !right.contains(&b'N'))
            .then(|| (k, left.to_vec()))
    })
}

/// Modal TSD length across copies with pinned ends, compared with the same
/// test at ends shifted into the flanks.
pub fn call_tsd(copies: &[&OrientedCopy], consensus: &[u8], p: &CopyParams) -> Option<Tsd> {
    let reach = p.tsd_max + p.tsd_null_shift;
    let ends = copies
        .iter()
        .filter_map(|c| pin_ends(c, consensus).map(|e| (c, e)))
        .filter(|(c, (s, e))| *s >= reach && e + reach <= c.seq.len())
        .collect::<Vec<_>>();
    if ends.len() < p.tsd_min_copies {
        return None;
    }
    let mut observed: Vec<Option<(usize, Vec<u8>)>> = Vec::new();
    let mut null: Vec<Option<usize>> = Vec::new();
    for (c, (s, e)) in &ends {
        observed.push(tsd_at(&c.seq, *s, *e, p));
        null.push(
            tsd_at(&c.seq, s - p.tsd_null_shift, e + p.tsd_null_shift, p).map(|(k, _)| k),
        );
    }
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for (k, _) in observed.iter().flatten() {
        *counts.entry(*k).or_default() += 1;
    }
    // the modal exact length (chance extensions of a TSD by a base or two are
    // common, so ties go to the shorter length); support tolerates ±1 bp
    let score = |k: usize| {
        (k.saturating_sub(1)..=k + 1)
            .map(|j| counts.get(&j).copied().unwrap_or(0))
            .sum::<usize>()
    };
    let mode = counts
        .keys()
        .copied()
        .max_by_key(|&k| (counts[&k], std::cmp::Reverse(k)))?;
    let n = ends.len() as f64;
    let support = score(mode) as f64 / n;
    let null_support = null
        .iter()
        .flatten()
        .filter(|&&k| k.abs_diff(mode) <= 1)
        .count() as f64
        / n;
    let mut motifs: HashMap<&[u8], usize> = HashMap::new();
    for (k, m) in observed.iter().flatten() {
        if *k == mode {
            *motifs.entry(m.as_slice()).or_default() += 1;
        }
    }
    let motif = motifs
        .into_iter()
        .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(a.0)))
        .map(|(m, _)| String::from_utf8_lossy(m).into_owned())
        .unwrap_or_default();
    (support >= 0.4 && support >= 2.0 * null_support).then_some(Tsd {
        length: mode,
        support,
        null_support,
        motif,
        copies_evaluated: ends.len(),
    })
}

/// Helitron insertion signature on genomic copies: the element starts with
/// TC and ends with CTRR, and sits between a host A and a host T without a
/// TSD (on the other strand: A|YYAG ... GA|T). Consensus ends are often a few
/// bases off, so each end is searched within a window and the best offset
/// kept; the same search far into the flanks gives the chance expectation.
#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct HelitronEnds {
    /// Fraction of copies showing the signature at both ends.
    pub support: f64,
    /// The same, searched 30 bp into the flanks.
    pub null_support: f64,
    pub copies: usize,
}

/// Bases searched on either side of each consensus end.
const END_WINDOW: isize = 8;
/// Shift into the flanks for the null search.
const END_NULL_SHIFT: isize = 30;

fn matches(seq: &[u8], at: isize, pattern: &[u8]) -> bool {
    let Ok(at) = usize::try_from(at) else {
        return false;
    };
    seq.get(at..at + pattern.len()).is_some_and(|s| {
        s.iter().zip(pattern).all(|(b, p)| match p {
            b'R' => matches!(b, b'A' | b'G'),
            b'Y' => matches!(b, b'C' | b'T'),
            _ => b == p,
        })
    })
}

/// Fraction of copies with `five` at the best 5' offset and `three` (ending
/// at the boundary) at the best 3' offset, with ends moved `shift` outward.
fn end_support(copies: &[&OrientedCopy], five: &[u8], three: &[u8], shift: isize) -> f64 {
    let at5 = |c: &OrientedCopy, k: isize| c.left as isize - shift + k;
    let at3 = |c: &OrientedCopy, k: isize| c.core_end() as isize + shift + k;
    let best = |hit: &dyn Fn(&OrientedCopy, isize) -> bool| {
        (-END_WINDOW..=END_WINDOW)
            .max_by_key(|&k| (copies.iter().filter(|c| hit(c, k)).count(), -k.abs()))
            .expect("non-empty window")
    };
    // `five` starts one host base before the element; `three` ends one after
    let hit5 = |c: &OrientedCopy, k: isize| matches(&c.seq, at5(c, k) - 1, five);
    let hit3 =
        |c: &OrientedCopy, k: isize| matches(&c.seq, at3(c, k) - three.len() as isize + 1, three);
    let (k5, k3) = (best(&hit5), best(&hit3));
    copies.iter().filter(|c| hit5(c, k5) && hit3(c, k3)).count() as f64 / copies.len() as f64
}

/// Helitron end signature on full-length copies (≥5 needed).
pub fn helitron_ends(copies: &[&OrientedCopy]) -> Option<HelitronEnds> {
    let reach = (END_NULL_SHIFT + END_WINDOW + 6) as usize;
    let usable = copies
        .iter()
        .copied()
        .filter(|c| c.left >= reach && c.right >= reach && c.core_end() - c.left >= 24)
        .collect::<Vec<_>>();
    if usable.len() < 5 {
        return None;
    }
    // A|TC ... CTRR|T, or on the other strand A|YYAG ... GA|T
    let support = |shift| {
        end_support(&usable, b"ATC", b"CTRRT", shift)
            .max(end_support(&usable, b"AYYAG", b"GAT", shift))
    };
    Some(HelitronEnds {
        support: support(0),
        null_support: support(END_NULL_SHIFT),
        copies: usable.len(),
    })
}

/// Mean pairwise identity of the outer flank windows (both sides) of copies.
pub fn flank_similarity(copies: &[OrientedCopy], window: usize) -> Option<f64> {
    let flanks = copies
        .iter()
        .filter(|c| c.left >= window && c.right >= window)
        .map(|c| {
            (
                &c.seq[c.left - window..c.left],
                &c.seq[c.core_end()..c.core_end() + window],
            )
        })
        .collect::<Vec<_>>();
    if flanks.len() < 2 {
        return None;
    }
    let identity = |a: &[u8], b: &[u8]| {
        a.iter()
            .zip(b)
            .filter(|(x, y)| x == y && **x != b'N')
            .count() as f64
            / a.len() as f64
    };
    let (mut total, mut pairs) = (0.0, 0usize);
    for i in 0..flanks.len() {
        for j in i + 1..flanks.len() {
            total += identity(flanks[i].0, flanks[j].0).max(identity(flanks[i].1, flanks[j].1));
            pairs += 1;
        }
    }
    Some(total / pairs as f64)
}

/// Expected TSD lengths for a RepeatMasker label; `Some(0..=0)` means no TSD.
pub fn expected_tsd(label: &str) -> Option<std::ops::RangeInclusive<usize>> {
    let checks: &[(&str, std::ops::RangeInclusive<usize>)] = &[
        ("DNA/TcMar", 2..=2),
        ("DNA/hAT", 5..=8),
        ("DNA/MULE", 8..=11),
        ("DNA/PIF-Harbinger", 3..=3),
        ("DNA/CMC", 2..=3),
        ("DNA/PiggyBac", 4..=4),
        ("DNA/P", 7..=8),
        ("DNA/Merlin", 8..=9),
        ("DNA/Maverick", 5..=6),
        ("DNA/Crypton", 0..=0),
        ("RC/Helitron", 0..=0),
        ("LTR/DIRS", 0..=0),
        ("LTR/Ngaro", 0..=0),
        ("LTR", 4..=6),
    ];
    checks
        .iter()
        .find(|(prefix, _)| label.starts_with(prefix))
        .map(|(_, range)| range.clone())
}

/// Superfamilies within `order` (e.g. `DNA`) whose expected TSD includes `length`.
pub fn superfamilies_with_tsd(order: &str, length: usize) -> Vec<&'static str> {
    [
        "DNA/TcMar",
        "DNA/hAT",
        "DNA/MULE",
        "DNA/PIF-Harbinger",
        "DNA/CMC",
        "DNA/PiggyBac",
        "DNA/P",
        "DNA/Merlin",
        "DNA/Maverick",
        "LTR",
    ]
    .into_iter()
    .filter(|sf| sf.starts_with(order) && expected_tsd(sf).is_some_and(|r| r.contains(&length)))
    .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn random(n: usize, seed: u64) -> Vec<u8> {
        // scramble the seed and discard early outputs: xorshift output from
        // small seeds is low-entropy at first, which creates spurious repeats
        let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        for _ in 0..16 {
            next();
        }
        (0..n)
            .map(|_| b"ACGT"[(next() >> 33) as usize % 4])
            .collect()
    }

    /// flank + TSD + element + TSD + flank, with the boundary shifted by `jitter`.
    fn copy(seed: u64, tsd_len: usize, jitter: usize) -> OrientedCopy {
        let left_flank = random(200, seed);
        let tsd = random(tsd_len, seed + 1000);
        let element = random(500, 42);
        let right_flank = random(200, seed + 2000);
        let seq = [left_flank, tsd.clone(), element, tsd, right_flank].concat();
        OrientedCopy {
            locus_id: format!("locus-{seed}"),
            left: 200 + tsd_len + jitter,
            right: 200 + tsd_len,
            seq,
        }
    }

    #[test]
    fn calls_modal_tsd_with_boundary_jitter() {
        let copies = (0..8)
            .map(|i| copy(i * 7 + 1, 8, (i % 3) as usize))
            .collect::<Vec<_>>();
        let refs = copies.iter().collect::<Vec<_>>();
        let tsd = call_tsd(&refs, &random(500, 42), &CopyParams::default()).unwrap();
        assert_eq!(tsd.length, 8);
        assert!(tsd.support >= 0.75, "{tsd:?}");
        assert!(tsd.null_support < 0.3);
        assert_eq!(tsd.copies_evaluated, 8);
    }

    #[test]
    fn pinned_ends_find_tsds_beyond_hits_that_stop_inside_the_element() {
        // BLAST hits stopping 9 bp inside each end, and a TSD with a point mutation
        let element = random(500, 42);
        let copies = (0..6u64)
            .map(|i| {
                let tsd = random(8, i + 1000);
                let mut right_tsd = tsd.clone();
                right_tsd[3] = if right_tsd[3] == b'A' { b'C' } else { b'A' };
                let seq = [
                    random(200, i + 1),
                    tsd,
                    element.clone(),
                    right_tsd,
                    random(200, i + 2000),
                ]
                .concat();
                OrientedCopy {
                    locus_id: i.to_string(),
                    left: 208 + 9,
                    right: 208 + 9,
                    seq,
                }
            })
            .collect::<Vec<_>>();
        let refs = copies.iter().collect::<Vec<_>>();
        assert_eq!(pin_ends(refs[0], &element), Some((208, 708)));
        let tsd = call_tsd(&refs, &element, &CopyParams::default()).unwrap();
        assert_eq!((tsd.length, tsd.copies_evaluated), (8, 6));
    }

    #[test]
    fn no_tsd_without_duplications() {
        let copies = (0..8).map(|i| copy(i * 7 + 1, 0, 0)).collect::<Vec<_>>();
        let refs = copies.iter().collect::<Vec<_>>();
        assert_eq!(call_tsd(&refs, &random(500, 42), &CopyParams::default()), None);
    }

    #[test]
    fn too_few_copies_give_no_call() {
        let copies = (0..2).map(|i| copy(i + 1, 8, 0)).collect::<Vec<_>>();
        assert_eq!(
            call_tsd(&copies.iter().collect::<Vec<_>>(), &random(500, 42), &CopyParams::default()),
            None
        );
    }

    #[test]
    fn flank_similarity_distinguishes_unique_and_shared_flanks() {
        let unique = (0..5).map(|i| copy(i * 11 + 3, 5, 0)).collect::<Vec<_>>();
        assert!(flank_similarity(&unique, 50).unwrap() < 0.4);
        // segmental-duplication-like: identical flanks around every copy
        let shared = (0..5)
            .map(|i| OrientedCopy {
                locus_id: format!("l{i}"),
                ..copy(99, 5, 0)
            })
            .collect::<Vec<_>>();
        assert!(flank_similarity(&shared, 50).unwrap() > 0.99);
    }

    #[test]
    fn expected_tsd_table() {
        assert_eq!(expected_tsd("DNA/TcMar-Mariner"), Some(2..=2));
        assert_eq!(expected_tsd("LTR/Gypsy"), Some(4..=6));
        assert_eq!(expected_tsd("RC/Helitron"), Some(0..=0));
        assert_eq!(expected_tsd("LINE/L1"), None);
        assert_eq!(
            superfamilies_with_tsd("DNA", 8),
            vec!["DNA/hAT", "DNA/MULE", "DNA/P", "DNA/Merlin"]
        );
    }

    #[test]
    fn reads_copy_recovery_family_directory() {
        let dir = tempfile::tempdir().unwrap();
        let d = dir.path();
        fs::write(
            d.join("loci.tsv"),
            "family_id\tlocus_id\tsubject_id\tstart\tend\textended_start\textended_end\tstrand\tn_fragments\tbest_bitscore\tcombined_query_coverage\tpercent_identity\trepresentative_score\tselected\texclusion_reason\n\
             f\tlocus-1\tc\t0\t10\t0\t10\t+\t1\t100\t1.0\t90\t1\ttrue\t\n\
             f\tlocus-2\tc\t0\t10\t0\t10\t-\t1\t100\t0.5\t80\t1\ttrue\t\n\
             f\tlocus-3\tc\t0\t10\t0\t10\t+\t1\t100\t1.0\t70\t1\tfalse\tduplicate_of:locus-1\n",
        )
        .unwrap();
        // locus-2 is reverse-complemented: genomic left 10 / right 20 -> oriented 20 / 10
        fs::write(
            d.join("copies.tsv"),
            "family_id\tlocus_id\tsubject_id\tgenomic_start\tgenomic_end\tstrand\trequested_extended_start\trequested_extended_end\tactual_extended_start\tactual_extended_end\tleft_flank_actual\tright_flank_actual\treverse_complemented\tselected\textraction_status\talignment_status\texclusion_reason\n\
             f\tlocus-1\tc\t0\t0\t+\t0\t0\t0\t0\t20\t20\tfalse\ttrue\tsuccess\tpending\t\n\
             f\tlocus-2\tc\t0\t0\t-\t0\t0\t0\t0\t10\t20\ttrue\ttrue\tsuccess\tpending\t\n",
        )
        .unwrap();
        fs::write(
            d.join("copies.extended.fa"),
            format!(
                ">locus_id=locus-1 family_id=f\n{}\n>locus_id=locus-2 family_id=f\n{}\n",
                "A".repeat(60),
                "C".repeat(60)
            ),
        )
        .unwrap();
        let e = read_copy_evidence(d, b"ACGT", &CopyParams::default())
            .unwrap()
            .unwrap();
        assert_eq!((e.loci, e.full_length), (2, 1));
        assert!((e.mean_identity - 85.0).abs() < 1e-9);
        let copies = oriented_copies(d).unwrap();
        assert_eq!((copies[1].left, copies[1].right), (20, 10));
        assert!(
            read_copy_evidence(&d.join("missing"), b"ACGT", &CopyParams::default())
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn helitron_ends_need_termini_and_host_a_t() {
        let mut x = 99u64;
        let mut rnd = |n: usize| {
            (0..n)
                .map(|_| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    b"ACGT"[(x >> 33) as usize % 4]
                })
                .collect::<Vec<u8>>()
        };
        let body = rnd(80);
        let mut make = |element: &[u8], host: (u8, u8), reverse: bool| {
            (0..8)
                .map(|i| {
                    let (mut left, mut right) = (rnd(70), rnd(70));
                    left[69] = host.0;
                    right[0] = host.1;
                    let mut seq = [left, element.to_vec(), right].concat();
                    if reverse {
                        seq = bio::alphabets::dna::revcomp(&seq);
                    }
                    OrientedCopy {
                        locus_id: format!("c{i}"),
                        seq,
                        // the consensus misses the first two and last three bases
                        left: 72,
                        right: 73,
                    }
                })
                .collect::<Vec<_>>()
        };
        let helitron = [b"TCTA".as_slice(), &body, b"GGCTAG"].concat();
        for reverse in [false, true] {
            let copies = make(&helitron, (b'A', b'T'), reverse);
            let e = helitron_ends(&copies.iter().collect::<Vec<_>>()).unwrap();
            assert_eq!(e.support, 1.0, "reverse={reverse}");
            assert!(e.null_support < 0.3);
        }
        // the same element between other host bases, or a DNA transposon end
        let copies = make(&helitron, (b'G', b'C'), false);
        assert!(
            helitron_ends(&copies.iter().collect::<Vec<_>>())
                .unwrap()
                .support
                < 0.3
        );
        let tir = [b"CAGTGG".as_slice(), &body, b"CCACTG"].concat();
        let copies = make(&tir, (b'A', b'T'), false);
        assert!(
            helitron_ends(&copies.iter().collect::<Vec<_>>())
                .unwrap()
                .support
                < 0.3
        );
    }
}
