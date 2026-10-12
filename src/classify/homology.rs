//! Protein homology of nucleotide consensi against a TE protein library.
//!
//! The search reports query (nucleotide) coordinates as 1-based inclusive, as
//! BLASTX and DIAMOND do; they are normalised so that start <= end.

use std::{collections::BTreeMap, path::Path, process::Command};

use csv::ReaderBuilder;
use serde::{Deserialize, Serialize};

use super::taxonomy::{class_from_protein_name, top_level};
use crate::{Error, ErrorKind, Result};

/// Columns requested from both BLASTX and DIAMOND.
pub const OUTFMT_FIELDS: &str =
    "qseqid sseqid qlen slen pident length qstart qend sstart send evalue bitscore";

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum HomologyBackend {
    Blastx,
    Diamond,
}

impl std::str::FromStr for HomologyBackend {
    type Err = String;
    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        match value {
            "blastx" => Ok(Self::Blastx),
            "diamond" => Ok(Self::Diamond),
            _ => Err(format!("unsupported homology backend: {value}")),
        }
    }
}

pub fn makeblastdb_protein_command(exe: impl AsRef<Path>, fasta: &Path, out: &Path) -> Command {
    let mut c = Command::new(exe.as_ref());
    c.arg("-in")
        .arg(fasta)
        .args(["-dbtype", "prot", "-out"])
        .arg(out);
    c
}

pub fn makeblastdb_nucleotide_command(exe: impl AsRef<Path>, fasta: &Path, out: &Path) -> Command {
    let mut c = Command::new(exe.as_ref());
    c.arg("-in")
        .arg(fasta)
        .args(["-dbtype", "nucl", "-out"])
        .arg(out);
    c
}

/// Consensus-versus-library nucleotide search. `dc-megablast` tolerates the
/// divergence between related families while staying fast on large libraries.
pub fn blastn_library_command(
    exe: impl AsRef<Path>,
    db: &Path,
    query: &Path,
    out: &Path,
    threads: usize,
) -> Command {
    let mut c = Command::new(exe.as_ref());
    c.args(["-task", "dc-megablast", "-db"])
        .arg(db)
        .arg("-query")
        .arg(query)
        .args([
            "-outfmt",
            &format!("6 {OUTFMT_FIELDS}"),
            "-evalue",
            "1e-10",
            "-max_target_seqs",
            "50",
            "-num_threads",
        ])
        .arg(threads.to_string())
        .arg("-out")
        .arg(out);
    c
}

pub fn blastx_command(
    exe: impl AsRef<Path>,
    db: &Path,
    query: &Path,
    out: &Path,
    threads: usize,
    evalue: f64,
) -> Command {
    let mut c = Command::new(exe.as_ref());
    c.arg("-db")
        .arg(db)
        .arg("-query")
        .arg(query)
        .args(["-outfmt", &format!("6 {OUTFMT_FIELDS}"), "-evalue"])
        .arg(evalue.to_string())
        .args(["-max_target_seqs", "50", "-seg", "yes", "-num_threads"])
        .arg(threads.to_string())
        .arg("-out")
        .arg(out);
    c
}

pub fn diamond_makedb_command(exe: impl AsRef<Path>, fasta: &Path, out: &Path) -> Command {
    let mut c = Command::new(exe.as_ref());
    c.args(["makedb", "--in"]).arg(fasta).arg("--db").arg(out);
    c
}

pub fn diamond_blastx_command(
    exe: impl AsRef<Path>,
    db: &Path,
    query: &Path,
    out: &Path,
    threads: usize,
    evalue: f64,
) -> Command {
    let mut c = Command::new(exe.as_ref());
    c.args(["blastx", "--more-sensitive", "--db"])
        .arg(db)
        .arg("--query")
        .arg(query)
        .arg("--out")
        .arg(out)
        .arg("--outfmt")
        .arg("6")
        .args(OUTFMT_FIELDS.split(' '))
        .arg("--evalue")
        .arg(evalue.to_string())
        .args(["--max-target-seqs", "50", "--threads"])
        .arg(threads.to_string());
    c
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ProteinHit {
    pub query_id: String,
    pub subject_id: String,
    /// Class inferred from the subject name, if any rule matched.
    pub class: Option<String>,
    pub query_length: u64,
    pub subject_length: u64,
    pub percent_identity: f64,
    pub query_start: u64,
    pub query_end: u64,
    /// `+` when the hit translates the consensus forward strand.
    pub strand: char,
    pub subject_start: u64,
    pub subject_end: u64,
    pub evalue: f64,
    pub bitscore: f64,
}

impl ProteinHit {
    fn from_row(row: &csv::StringRecord) -> Result<Self> {
        if row.len() != 12 {
            return Err(Error::new(ErrorKind::GenericCli(format!(
                "Expected 12 protein-search columns ({OUTFMT_FIELDS}), found {} in {:?}",
                row.len(),
                row
            ))));
        }
        let qstart: u64 = row[6].parse()?;
        let qend: u64 = row[7].parse()?;
        Ok(Self {
            query_id: row[0].to_owned(),
            subject_id: row[1].to_owned(),
            class: class_from_protein_name(&row[1]),
            query_length: row[2].parse()?,
            subject_length: row[3].parse()?,
            percent_identity: row[4].parse()?,
            query_start: qstart.min(qend),
            query_end: qstart.max(qend),
            strand: if qend >= qstart { '+' } else { '-' },
            subject_start: row[8].parse()?,
            subject_end: row[9].parse()?,
            evalue: row[10].parse()?,
            bitscore: row[11].parse()?,
        })
    }
}

pub fn parse_protein_hits(path: &Path) -> Result<Vec<ProteinHit>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .comment(Some(b'#'))
        .from_path(path)?;
    reader
        .records()
        .map(|row| ProteinHit::from_row(&row?))
        .collect()
}

/// Aggregated support for one class from all hits of one family.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ClassSupport {
    pub class: String,
    pub best_bitscore: f64,
    pub best_subject: String,
    pub best_identity: f64,
    /// Fraction of the best subject covered by its hit.
    pub best_subject_coverage: f64,
    /// Fraction of the consensus covered by hits to this class.
    pub query_coverage: f64,
    pub n_hits: usize,
    /// Merged consensus intervals (0-based half-open) covered by this class.
    pub intervals: Vec<(u64, u64)>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct HomologySummary {
    /// Class supports in descending best-bitscore order.
    pub classes: Vec<ClassSupport>,
    /// Hits whose subject name matched no taxonomy rule.
    pub unclassified_hits: usize,
    pub best_unclassified_subject: Option<String>,
    pub best_unclassified_bitscore: f64,
    /// Top-level classes supported in largely separate consensus regions.
    pub chimeric_classes: Vec<String>,
}

impl HomologySummary {
    pub fn best(&self) -> Option<&ClassSupport> {
        self.classes.first()
    }
}

/// Summarise the hits of a single family. Hits below `min_bitscore` are ignored.
pub fn summarise(hits: &[&ProteinHit], min_bitscore: f64) -> HomologySummary {
    let mut by_class: BTreeMap<String, Vec<&ProteinHit>> = BTreeMap::new();
    let mut summary = HomologySummary::default();
    for &hit in hits.iter().filter(|h| h.bitscore >= min_bitscore) {
        match &hit.class {
            Some(class) => by_class.entry(class.clone()).or_default().push(hit),
            None => {
                summary.unclassified_hits += 1;
                if hit.bitscore > summary.best_unclassified_bitscore {
                    summary.best_unclassified_bitscore = hit.bitscore;
                    summary.best_unclassified_subject = Some(hit.subject_id.clone());
                }
            }
        }
    }
    for (class, class_hits) in by_class {
        let best = class_hits
            .iter()
            .max_by(|a, b| {
                a.bitscore
                    .total_cmp(&b.bitscore)
                    .then_with(|| b.subject_id.cmp(&a.subject_id))
            })
            .expect("non-empty group");
        let intervals = merge_intervals(
            class_hits
                .iter()
                .map(|h| (h.query_start - 1, h.query_end))
                .collect(),
        );
        let covered: u64 = intervals.iter().map(|(s, e)| e - s).sum();
        summary.classes.push(ClassSupport {
            class,
            best_bitscore: best.bitscore,
            best_subject: best.subject_id.clone(),
            best_identity: best.percent_identity,
            best_subject_coverage: (best.subject_start.abs_diff(best.subject_end) + 1) as f64
                / best.subject_length.max(1) as f64,
            query_coverage: covered as f64 / best.query_length.max(1) as f64,
            n_hits: class_hits.len(),
            intervals,
        });
    }
    summary.classes.sort_by(|a, b| {
        b.best_bitscore
            .total_cmp(&a.best_bitscore)
            .then_with(|| a.class.cmp(&b.class))
    });
    summary.chimeric_classes = chimeric_classes(&summary.classes);
    summary
}

type Interval = (u64, u64);

/// Top-level classes (LTR, LINE, DNA, ...) that each have strong support in
/// mostly non-overlapping regions of the consensus. Superfamilies of the same
/// class are pooled first, because related proteins routinely cross-hit.
fn chimeric_classes(classes: &[ClassSupport]) -> Vec<String> {
    let mut pooled: BTreeMap<&str, (f64, Vec<Interval>)> = BTreeMap::new();
    for support in classes {
        let entry = pooled.entry(top_level(&support.class)).or_default();
        entry.0 = entry.0.max(support.best_bitscore);
        entry.1.extend(support.intervals.iter().copied());
    }
    let mut pooled = pooled
        .into_iter()
        .map(|(class, (score, intervals))| (class, score, merge_intervals(intervals)))
        .collect::<Vec<_>>();
    pooled.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(b.0)));
    let Some((first, first_score, first_iv)) = pooled.first().cloned() else {
        return Vec::new();
    };
    let mut result = vec![first.to_owned()];
    for (class, score, iv) in pooled.iter().skip(1) {
        let own: u64 = iv.iter().map(|(s, e)| e - s).sum();
        let shared = overlap(&first_iv, iv);
        if *score >= 0.5 * first_score && own > 0 && (shared as f64) < 0.5 * own as f64 {
            result.push((*class).to_owned());
        }
    }
    if result.len() > 1 {
        result
    } else {
        Vec::new()
    }
}

fn merge_intervals(mut intervals: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    intervals.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::new();
    for (s, e) in intervals {
        match merged.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => merged.push((s, e)),
        }
    }
    merged
}

fn overlap(a: &[(u64, u64)], b: &[(u64, u64)]) -> u64 {
    let mut total = 0;
    for &(s1, e1) in a {
        for &(s2, e2) in b {
            total += e1.min(e2).saturating_sub(s1.max(s2));
        }
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(subject: &str, qs: u64, qe: u64, bits: f64) -> ProteinHit {
        ProteinHit::from_row(&csv::StringRecord::from(vec![
            "fam".to_string(),
            subject.to_string(),
            "3000".to_string(),
            "500".to_string(),
            "45.0".to_string(),
            "300".to_string(),
            qs.to_string(),
            qe.to_string(),
            "1".to_string(),
            "300".to_string(),
            "1e-30".to_string(),
            bits.to_string(),
        ]))
        .unwrap()
    }

    #[test]
    fn parses_reverse_frame_hits() {
        let h = hit("Gypsy-1_xx", 900, 1, 200.0);
        assert_eq!((h.query_start, h.query_end, h.strand), (1, 900, '-'));
        assert_eq!(h.class.as_deref(), Some("LTR/Gypsy"));
    }

    #[test]
    fn summarises_classes_and_coverage() {
        let hits = [
            hit("Gypsy-1_xx", 1, 900, 300.0),
            hit("Gypsy-2_xx", 801, 1500, 150.0),
            hit("Copia-1_xx", 1, 600, 90.0),
            hit("Itch-1_xx", 1, 300, 60.0),
            hit("Gypsy-3_xx", 2000, 2100, 10.0),
        ];
        let refs = hits.iter().collect::<Vec<_>>();
        let s = summarise(&refs, 30.0);
        let best = s.best().unwrap();
        assert_eq!(best.class, "LTR/Gypsy");
        assert_eq!(best.intervals, vec![(0, 1500)]);
        assert!((best.query_coverage - 0.5).abs() < 1e-9);
        assert_eq!(s.unclassified_hits, 1);
        assert!(
            s.chimeric_classes.is_empty(),
            "same top-level class is not chimeric"
        );
    }

    #[test]
    fn detects_cross_class_chimeras() {
        let hits = [
            hit("Gypsy-1_xx", 1, 900, 300.0),
            hit("hAT-1_xx", 2001, 2900, 250.0),
        ];
        let refs = hits.iter().collect::<Vec<_>>();
        assert_eq!(summarise(&refs, 30.0).chimeric_classes, vec!["LTR", "DNA"]);
    }

    #[test]
    fn command_lines_request_the_parsed_columns() {
        let c = blastx_command(
            "blastx",
            Path::new("db"),
            Path::new("q"),
            Path::new("o"),
            2,
            1e-5,
        );
        let args = c
            .get_args()
            .map(|a| a.to_string_lossy().into_owned())
            .collect::<Vec<_>>();
        assert!(args.contains(&format!("6 {OUTFMT_FIELDS}")));
        let c = diamond_blastx_command(
            "diamond",
            Path::new("db"),
            Path::new("q"),
            Path::new("o"),
            2,
            1e-5,
        );
        let args = c
            .get_args()
            .map(|a| a.to_string_lossy().into_owned())
            .collect::<Vec<_>>();
        let i = args.iter().position(|a| a == "--outfmt").unwrap();
        assert_eq!(args[i + 1..i + 14].join(" "), format!("6 {OUTFMT_FIELDS}"));
    }
}
