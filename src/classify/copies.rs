//! Evidence from recovered genomic copies (`rep recover-copies` output).
//!
//! For each family: copy number and completeness from `loci.tsv`, divergence
//! from the consensus, target-site duplications (TSDs) at the boundaries of
//! full-length copies, and whether the flanks of different copies resemble
//! each other (expected for segmental duplications or truncated consensi, not
//! for transposable elements).

use std::{collections::HashMap, path::Path};

use bio::io::fasta;
use serde::Serialize;

use crate::Result;

#[derive(Clone, Debug, Serialize)]
pub struct CopyParams {
    /// Combined consensus coverage for a copy to count as full-length.
    pub full_length_coverage: f64,
    /// Bases of boundary uncertainty searched on each side for a TSD.
    pub tsd_slack: usize,
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
            tsd_slack: 3,
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
}

/// Read copy evidence for one family from its copy-recovery directory.
/// Returns `None` when the directory does not exist.
pub fn read_copy_evidence(dir: &Path, params: &CopyParams) -> Result<Option<CopyEvidence>> {
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
        if reason.starts_with("duplicate_of") || reason == "ambiguous_strand" {
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
    evidence.tsd = call_tsd(&complete, params);
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

/// Longest exact duplication flanking the boundaries `start` and `end`,
/// allowing each boundary to move by up to `slack` bases.
fn best_tsd(seq: &[u8], start: usize, end: usize, p: &CopyParams) -> Option<(usize, Vec<u8>)> {
    for k in (p.tsd_min..=p.tsd_max).rev() {
        for dl in -(p.tsd_slack as isize)..=p.tsd_slack as isize {
            for dr in -(p.tsd_slack as isize)..=p.tsd_slack as isize {
                let s = start as isize + dl;
                let e = end as isize + dr;
                if s < k as isize || e < 0 || e as usize + k > seq.len() || e <= s {
                    continue;
                }
                let (s, e) = (s as usize, e as usize);
                let left = &seq[s - k..s];
                if left == &seq[e..e + k] && !left.contains(&b'N') {
                    return Some((k, left.to_vec()));
                }
            }
        }
    }
    None
}

/// Modal TSD length across copies, compared with the same search at
/// boundaries shifted into the flanks.
pub fn call_tsd(copies: &[&OrientedCopy], p: &CopyParams) -> Option<Tsd> {
    let usable = copies
        .iter()
        .filter(|c| {
            c.left >= p.tsd_max + p.tsd_slack + p.tsd_null_shift
                && c.right >= p.tsd_max + p.tsd_slack + p.tsd_null_shift
        })
        .collect::<Vec<_>>();
    if usable.len() < p.tsd_min_copies {
        return None;
    }
    let mut observed: Vec<Option<(usize, Vec<u8>)>> = Vec::new();
    let mut null: Vec<Option<usize>> = Vec::new();
    for c in &usable {
        observed.push(best_tsd(&c.seq, c.left, c.core_end(), p));
        null.push(
            best_tsd(
                &c.seq,
                c.left - p.tsd_null_shift,
                c.core_end() + p.tsd_null_shift,
                p,
            )
            .map(|(k, _)| k),
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
    let n = usable.len() as f64;
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
        copies_evaluated: usable.len(),
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
        let tsd = call_tsd(&refs, &CopyParams::default()).unwrap();
        assert_eq!(tsd.length, 8);
        assert!(tsd.support >= 0.75, "{tsd:?}");
        assert!(tsd.null_support < 0.3);
        assert_eq!(tsd.copies_evaluated, 8);
    }

    #[test]
    fn no_tsd_without_duplications() {
        let copies = (0..8).map(|i| copy(i * 7 + 1, 0, 0)).collect::<Vec<_>>();
        let refs = copies.iter().collect::<Vec<_>>();
        assert_eq!(call_tsd(&refs, &CopyParams::default()), None);
    }

    #[test]
    fn too_few_copies_give_no_call() {
        let copies = (0..2).map(|i| copy(i + 1, 8, 0)).collect::<Vec<_>>();
        assert_eq!(
            call_tsd(&copies.iter().collect::<Vec<_>>(), &CopyParams::default()),
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
        let e = read_copy_evidence(d, &CopyParams::default())
            .unwrap()
            .unwrap();
        assert_eq!((e.loci, e.full_length), (2, 1));
        assert!((e.mean_identity - 85.0).abs() < 1e-9);
        let copies = oriented_copies(d).unwrap();
        assert_eq!((copies[1].left, copies[1].right), (20, 10));
        assert!(
            read_copy_evidence(&d.join("missing"), &CopyParams::default())
                .unwrap()
                .is_none()
        );
    }
}
