use std::path::Path;

use csv::ReaderBuilder;
use serde::{Deserialize, Serialize};

use crate::{Error, ErrorKind, Result};

#[derive(Clone, Copy, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd, Serialize)]
pub enum Strand {
    Plus,
    Minus,
}

impl Strand {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Plus => "+",
            Self::Minus => "-",
        }
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct Hit {
    pub family_id: String,
    pub query_id: String,
    pub subject_id: String,
    pub query_start: u64,
    pub query_end: u64,
    pub query_length: u64,
    pub subject_start: u64,
    pub subject_end: u64,
    pub strand: Strand,
    pub alignment_length: u64,
    pub matches: u64,
    pub mismatches: u64,
    pub gap_opens: u64,
    pub percent_identity: f64,
    pub query_coverage: f64,
    pub bitscore: f64,
    pub evalue: f64,
}

impl Hit {
    fn from_row(row: &csv::StringRecord) -> Result<Self> {
        if row.len() != 12 && row.len() != 13 {
            return Err(Error::new(ErrorKind::GenericCli(format!(
                "Expected 12 or 13 copy-recovery BLAST columns, found {} in record {:?}",
                row.len(),
                row
            ))));
        }
        let qstart: u64 = row[7].parse()?;
        let qend: u64 = row[8].parse()?;
        let qlen: u64 = row[2].parse()?;
        let raw_start: u64 = row[9].parse()?;
        let raw_end: u64 = row[10].parse()?;
        let alignment_length: u64 = row[4].parse()?;
        let mismatches: u64 = row[5].parse()?;
        let gap_opens: u64 = row[6].parse()?;
        let percent_identity: f64 = row[3].parse()?;
        Ok(Self {
            family_id: row[0].to_owned(),
            query_id: row[0].to_owned(),
            subject_id: row[1].to_owned(),
            query_start: qstart.min(qend),
            query_end: qstart.max(qend),
            query_length: qlen,
            subject_start: raw_start.min(raw_end),
            subject_end: raw_start.max(raw_end),
            strand: if raw_end >= raw_start {
                Strand::Plus
            } else {
                Strand::Minus
            },
            alignment_length,
            matches: alignment_length.saturating_sub(mismatches + gap_opens),
            mismatches,
            gap_opens,
            percent_identity,
            query_coverage: if qlen == 0 {
                0.0
            } else {
                (qstart.abs_diff(qend) + 1) as f64 / qlen as f64
            },
            evalue: row[if row.len() == 12 { 10 } else { 11 }].parse()?,
            bitscore: row[if row.len() == 12 { 11 } else { 12 }].parse()?,
        })
    }
}

pub fn parse_blast_file(path: impl AsRef<Path>) -> Result<Vec<Hit>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .comment(Some(b'#'))
        .from_path(path)?;
    reader.records().map(|row| Hit::from_row(&row?)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_and_normalizes_minus_strand_row() {
        let row = csv::StringRecord::from(vec![
            "fam", "chr1", "100", "95.0", "50", "2", "1", "10", "59", "900", "851", "1e-20", "80",
        ]);
        let hit = Hit::from_row(&row).unwrap();
        assert_eq!(hit.strand, Strand::Minus);
        assert_eq!((hit.subject_start, hit.subject_end), (851, 900));
        assert_eq!(hit.matches, 47);
        assert_eq!(hit.query_coverage, 0.5);
    }
}
