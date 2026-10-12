//! Sequence-intrinsic structural evidence for a consensus, with no external tools.
//!
//! All coordinates are 0-based half-open on the forward strand of the
//! consensus. Input sequences are expected to be upper case.

use bio::alignment::pairwise::Aligner;
use bio::alignment::AlignmentOperation;
use bio::alphabets::dna;
use serde::Serialize;

/// Parameters for the structural detectors.
#[derive(Clone, Debug, Serialize)]
pub struct StructureParams {
    pub ltr_window: usize,
    pub ltr_min_length: usize,
    pub tir_window: usize,
    pub tir_min_length: usize,
    pub terminal_min_identity: f64,
    /// Bases tolerated between a terminal repeat and the consensus end.
    pub terminal_slack: usize,
    pub tandem_max_period: usize,
    pub tandem_sample: usize,
    pub tandem_min_identity: f64,
    pub min_orf_aa: usize,
}

impl Default for StructureParams {
    fn default() -> Self {
        Self {
            ltr_window: 2500,
            ltr_min_length: 80,
            tir_window: 500,
            tir_min_length: 10,
            terminal_min_identity: 0.8,
            terminal_slack: 10,
            tandem_max_period: 1000,
            tandem_sample: 3000,
            tandem_min_identity: 0.7,
            min_orf_aa: 100,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct TerminalRepeat {
    /// Length of the 5' copy of the repeat.
    pub length: usize,
    pub identity: f64,
    /// Distance of the 5' copy from the consensus start.
    pub five_prime_offset: usize,
    /// Distance of the 3' copy from the consensus end.
    pub three_prime_offset: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Tail {
    /// Poly-A at the 3' end (non-LTR retrotransposon / SINE orientation).
    PolyA3,
    /// Poly-T at the 5' end (the same element in reverse orientation).
    PolyT5,
}

impl Tail {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::PolyA3 => "polyA_3prime",
            Self::PolyT5 => "polyT_5prime",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Tandem {
    pub period: usize,
    pub identity: f64,
    pub copies: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Composition {
    pub gc: f64,
    pub max_base_fraction: f64,
    pub n_fraction: f64,
    pub low_complexity: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Orf {
    /// `+` or `-` relative to the consensus.
    pub strand: char,
    pub frame: u8,
    pub start: usize,
    pub end: usize,
    pub protein: String,
}

impl Orf {
    pub fn aa_len(&self) -> usize {
        self.protein.len()
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct StructuralEvidence {
    pub length: usize,
    pub ltr: Option<TerminalRepeat>,
    pub tir: Option<TerminalRepeat>,
    pub tail: Option<Tail>,
    pub tandem: Option<Tandem>,
    pub composition: Composition,
    pub orfs: Vec<Orf>,
    pub five_prime: String,
    pub three_prime: String,
    /// Starts with TG and ends with CA, as LTRs and LTR elements do.
    pub tg_ca_termini: bool,
    /// Starts with TC and ends with CTRR, as Helitrons do.
    pub helitron_termini: bool,
}

impl StructuralEvidence {
    pub fn longest_orf_aa(&self) -> usize {
        self.orfs.iter().map(Orf::aa_len).max().unwrap_or(0)
    }
}

pub fn analyse(seq: &[u8], params: &StructureParams) -> StructuralEvidence {
    let five_prime = String::from_utf8_lossy(&seq[..seq.len().min(4)]).into_owned();
    let three_prime = String::from_utf8_lossy(&seq[seq.len().saturating_sub(4)..]).into_owned();
    let helitron_termini = seq.starts_with(b"TC")
        && seq.len() >= 4
        && &seq[seq.len() - 4..seq.len() - 2] == b"CT"
        && seq[seq.len() - 2..]
            .iter()
            .all(|b| matches!(b, b'A' | b'G'));
    let tandem = tandem(
        seq,
        params.tandem_max_period,
        params.tandem_sample,
        params.tandem_min_identity,
    );
    StructuralEvidence {
        length: seq.len(),
        // any tandem array trivially has matching ends, so the periodicity
        // explains the direct repeat rather than an LTR
        ltr: tandem
            .is_none()
            .then(|| {
                terminal_repeat(
                    seq,
                    false,
                    params.ltr_window,
                    params.ltr_min_length,
                    params.terminal_min_identity,
                    params.terminal_slack,
                )
            })
            .flatten(),
        tir: terminal_repeat(
            seq,
            true,
            params.tir_window,
            params.tir_min_length,
            params.terminal_min_identity,
            params.terminal_slack,
        ),
        tail: tail(seq),
        tandem,
        composition: composition(seq),
        orfs: orfs(seq, params.min_orf_aa),
        tg_ca_termini: seq.starts_with(b"TG") && seq.ends_with(b"CA"),
        helitron_termini,
        five_prime,
        three_prime,
    }
}

/// Find a repeat anchored at both consensus ends, either direct (LTR-like)
/// or inverted (TIR-like), by local alignment of the terminal windows.
pub fn terminal_repeat(
    seq: &[u8],
    inverted: bool,
    window: usize,
    min_length: usize,
    min_identity: f64,
    slack: usize,
) -> Option<TerminalRepeat> {
    let n = seq.len();
    if n < 2 * min_length {
        return None;
    }
    let w = window.min(n / 2);
    let x = &seq[..w];
    let suffix = &seq[n - w..];
    let y = if inverted {
        dna::revcomp(suffix)
    } else {
        suffix.to_vec()
    };
    let score = |a: u8, b: u8| if a == b && a != b'N' { 1i32 } else { -1i32 };
    let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, score);
    let aln = aligner.local(x, &y);
    if aln.xstart > slack {
        return None;
    }
    // y is the suffix itself (direct) or its reverse complement (inverted):
    // in both cases the end of the consensus must be within `slack` of the
    // corresponding end of the aligned y interval.
    let three_prime_offset = if inverted { aln.ystart } else { w - aln.yend };
    if three_prime_offset > slack {
        return None;
    }
    let (mut matches, mut columns) = (0usize, 0usize);
    for op in &aln.operations {
        match op {
            AlignmentOperation::Match => {
                matches += 1;
                columns += 1;
            }
            AlignmentOperation::Subst | AlignmentOperation::Del | AlignmentOperation::Ins => {
                columns += 1
            }
            _ => {}
        }
    }
    let length = aln.xend - aln.xstart;
    let identity = if columns == 0 {
        0.0
    } else {
        matches as f64 / columns as f64
    };
    (length >= min_length && identity >= min_identity).then_some(TerminalRepeat {
        length,
        identity,
        five_prime_offset: aln.xstart,
        three_prime_offset,
    })
}

/// A run of at least 11 A in 12 bases within the last 60 bases (or the same
/// for T at the 5' end), ignoring terminal Ns.
pub fn tail(seq: &[u8]) -> Option<Tail> {
    fn rich(region: &[u8], base: u8) -> bool {
        region
            .windows(12)
            .any(|w| w.iter().filter(|&&b| b == base).count() >= 11)
    }
    let trimmed_end = seq.len() - seq.iter().rev().take_while(|&&b| b == b'N').count();
    let trimmed_start = seq.iter().take_while(|&&b| b == b'N').count();
    if trimmed_end <= trimmed_start {
        return None;
    }
    let core = &seq[trimmed_start..trimmed_end];
    if rich(&core[core.len().saturating_sub(60)..], b'A') {
        Some(Tail::PolyA3)
    } else if rich(&core[..core.len().min(60)], b'T') {
        Some(Tail::PolyT5)
    } else {
        None
    }
}

/// Detect tandem periodicity: the fraction of positions identical to the base
/// one period downstream. The smallest period within 0.02 of the best score is
/// reported so that multiples of the true period are not preferred.
pub fn tandem(seq: &[u8], max_period: usize, sample: usize, min_identity: f64) -> Option<Tandem> {
    let s = &seq[..seq.len().min(sample)];
    let max_p = max_period.min(s.len() / 3);
    if max_p == 0 {
        return None;
    }
    let scores = (1..=max_p)
        .map(|p| {
            let same = s
                .iter()
                .zip(&s[p..])
                .filter(|(a, b)| a == b && **a != b'N')
                .count();
            (p, same as f64 / (s.len() - p) as f64)
        })
        .collect::<Vec<_>>();
    let best = scores.iter().map(|&(_, x)| x).fold(0.0, f64::max);
    let (period, identity) = *scores.iter().find(|&&(_, x)| x >= best - 0.02)?;
    (identity >= min_identity).then(|| Tandem {
        period,
        identity,
        copies: seq.len() as f64 / period as f64,
    })
}

pub fn composition(seq: &[u8]) -> Composition {
    let mut counts = [0usize; 5];
    for &b in seq {
        counts[match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        }] += 1;
    }
    let called = counts[..4].iter().sum::<usize>().max(1) as f64;
    let f = counts[..4]
        .iter()
        .map(|&c| c as f64 / called)
        .collect::<Vec<_>>();
    let max_base_fraction = f.iter().copied().fold(0.0, f64::max);
    let max_pair = [(0, 3), (1, 2), (0, 2), (1, 3), (0, 1), (2, 3)]
        .iter()
        .map(|&(a, b)| f[a] + f[b])
        .fold(0.0, f64::max);
    Composition {
        gc: f[1] + f[2],
        max_base_fraction,
        n_fraction: counts[4] as f64 / seq.len().max(1) as f64,
        low_complexity: max_base_fraction >= 0.7 || max_pair >= 0.95,
    }
}

/// Stop-to-stop open reading frames of at least `min_aa` residues in all six
/// frames. ATG is not required because consensi of old families are often
/// degraded or truncated at the 5' end.
pub fn orfs(seq: &[u8], min_aa: usize) -> Vec<Orf> {
    let n = seq.len();
    let rc = dna::revcomp(seq);
    let mut out = Vec::new();
    for (strand, s) in [('+', seq), ('-', rc.as_slice())] {
        for frame in 0..3usize {
            let mut protein = String::new();
            let mut start = frame;
            let mut pos = frame;
            while pos + 3 <= n {
                let aa = translate(&s[pos..pos + 3]);
                if aa == '*' {
                    push_orf(&mut out, strand, frame, start, pos, &mut protein, n, min_aa);
                    start = pos + 3;
                } else {
                    protein.push(aa);
                }
                pos += 3;
            }
            push_orf(&mut out, strand, frame, start, pos, &mut protein, n, min_aa);
        }
    }
    out.sort_by(|a, b| b.aa_len().cmp(&a.aa_len()).then(a.start.cmp(&b.start)));
    out
}

#[allow(clippy::too_many_arguments)]
fn push_orf(
    out: &mut Vec<Orf>,
    strand: char,
    frame: usize,
    start: usize,
    end: usize,
    protein: &mut String,
    n: usize,
    min_aa: usize,
) {
    let p = std::mem::take(protein);
    if p.len() < min_aa {
        return;
    }
    // convert reverse-strand coordinates back to the forward strand
    let (start, end) = if strand == '+' {
        (start, end)
    } else {
        (n - end, n - start)
    };
    out.push(Orf {
        strand,
        frame: frame as u8,
        start,
        end,
        protein: p,
    });
}

/// Standard genetic code; codons containing anything other than ACGT give `X`.
pub fn translate(codon: &[u8]) -> char {
    const AA: &[u8; 64] = b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
    let index = |b: u8| match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    };
    match (index(codon[0]), index(codon[1]), index(codon[2])) {
        (Some(a), Some(b), Some(c)) => AA[a * 16 + b * 4 + c] as char,
        _ => 'X',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize]
            })
            .collect()
    }

    #[test]
    fn genetic_code() {
        assert_eq!(translate(b"ATG"), 'M');
        assert_eq!(translate(b"TAA"), '*');
        assert_eq!(translate(b"TGG"), 'W');
        assert_eq!(translate(b"GGN"), 'X');
        assert_eq!(translate(b"TTT"), 'F');
        assert_eq!(translate(b"GCA"), 'A');
    }

    #[test]
    fn detects_ltr_like_direct_repeats() {
        let ltr = random(300, 1);
        let seq = [ltr.clone(), random(2000, 2), ltr].concat();
        let tr = terminal_repeat(&seq, false, 2500, 80, 0.8, 10).unwrap();
        assert!(tr.length >= 295, "{tr:?}");
        assert!(terminal_repeat(&seq, true, 500, 10, 0.8, 10).is_none());
    }

    #[test]
    fn detects_terminal_inverted_repeats() {
        let tir = random(40, 3);
        let seq = [tir.clone(), random(1500, 4), dna::revcomp(&tir)].concat();
        let tr = terminal_repeat(&seq, true, 500, 10, 0.8, 10).unwrap();
        assert!(tr.length >= 38, "{tr:?}");
        assert!(terminal_repeat(&seq, false, 2500, 80, 0.8, 10).is_none());
    }

    #[test]
    fn random_sequence_has_no_structure() {
        let seq = random(3000, 5);
        let e = analyse(&seq, &StructureParams::default());
        assert!(e.ltr.is_none() && e.tir.is_none() && e.tail.is_none() && e.tandem.is_none());
        assert!(!e.composition.low_complexity);
    }

    #[test]
    fn detects_tails() {
        let body = random(500, 6);
        assert_eq!(
            tail(&[body.clone(), b"AAAAAAAAAAAAAAANN".to_vec()].concat()),
            Some(Tail::PolyA3)
        );
        assert_eq!(
            tail(&[b"TTTTTTTTTTTTTT".to_vec(), body.clone()].concat()),
            Some(Tail::PolyT5)
        );
        assert_eq!(tail(&body), None);
    }

    #[test]
    fn tandem_arrays_are_not_reported_as_ltrs() {
        let e = analyse(&random(178, 8).repeat(12), &StructureParams::default());
        assert!(e.tandem.is_some() && e.ltr.is_none());
    }

    #[test]
    fn detects_satellite_period() {
        let unit = random(171, 7);
        let seq = unit.repeat(10);
        let t = tandem(&seq, 1000, 3000, 0.7).unwrap();
        assert_eq!(t.period, 171);
        assert!(t.identity > 0.99);
        assert_eq!(
            tandem(b"CACACACACACACACACACA", 1000, 3000, 0.7)
                .unwrap()
                .period,
            2
        );
    }

    #[test]
    fn finds_orfs_on_both_strands_in_forward_coordinates() {
        // 120 codons without stops: GCT repeated, flanked by stops
        let coding = [b"TAA".to_vec(), b"GCT".repeat(120), b"TAA".to_vec()].concat();
        let plus = orfs(&coding, 100);
        assert!(plus
            .iter()
            .any(|o| o.strand == '+' && o.aa_len() == 120 && o.start == 3 && o.end == 363));
        let minus = orfs(&dna::revcomp(&coding), 100);
        assert!(minus
            .iter()
            .any(|o| o.strand == '-' && o.aa_len() == 120 && o.start == 3 && o.end == 363));
    }

    #[test]
    fn flags_low_complexity_and_termini() {
        assert!(composition(b"AAAAAAATAAAAAAAAAATA").low_complexity);
        let e = analyse(b"TCAGGCTAGCTTACGATCTAG", &StructureParams::default());
        assert!(e.helitron_termini);
        let e = analyse(b"TGAGGCTAGCTTACGATCTCA", &StructureParams::default());
        assert!(e.tg_ca_termini);
    }
}
