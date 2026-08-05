use crate::{Error, ErrorKind, Result};
use bio::alphabets::dna;
use bio::io::fasta;
use std::collections::BTreeMap;
use std::path::Path;

pub trait SequenceProvider {
    fn sequence(&self, contig: &str, start: u64, end: u64) -> Result<Vec<u8>>;
    fn contig_length(&self, contig: &str) -> Result<u64>;
}

pub struct FastaSequenceProvider {
    sequences: BTreeMap<String, Vec<u8>>,
}
impl FastaSequenceProvider {
    pub fn from_file(path: &Path) -> Result<Self> {
        let mut sequences = BTreeMap::new();
        let reader = fasta::Reader::from_file(path)?;
        for rec in reader.records() {
            let rec = rec?;
            if sequences
                .insert(rec.id().to_owned(), rec.seq().to_vec())
                .is_some()
            {
                return Err(Error::new(ErrorKind::GenericCli(format!(
                    "duplicate genome contig {}",
                    rec.id()
                ))));
            }
        }
        Ok(Self { sequences })
    }
}
impl SequenceProvider for FastaSequenceProvider {
    fn sequence(&self, contig: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        let seq = self.sequences.get(contig).ok_or_else(|| {
            Error::new(ErrorKind::GenericCli(format!(
                "unknown genome contig {contig}"
            )))
        })?;
        if start > end || end > seq.len() as u64 {
            return Err(Error::new(ErrorKind::GenericCli(format!(
                "invalid interval {contig}:{start}-{end}"
            ))));
        }
        Ok(seq[start as usize..end as usize].to_vec())
    }
    fn contig_length(&self, contig: &str) -> Result<u64> {
        self.sequences
            .get(contig)
            .map(|s| s.len() as u64)
            .ok_or_else(|| {
                Error::new(ErrorKind::GenericCli(format!(
                    "unknown genome contig {contig}"
                )))
            })
    }
}

pub fn clip_interval(start: u64, end: u64, length: usize, flank: usize) -> (usize, usize) {
    let s = start.saturating_sub(flank as u64) as usize;
    let e = end.saturating_add(flank as u64).min(length as u64) as usize;
    (s.min(e), e)
}
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    dna::revcomp(seq)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    #[test]
    fn clips_both_contig_boundaries() {
        assert_eq!(clip_interval(1, 4, 10, 5), (0, 9));
        assert_eq!(clip_interval(7, 10, 10, 5), (2, 10));
    }
    #[test]
    fn reverse_complements() {
        assert_eq!(reverse_complement(b"ACGTN"), b"NACGT".to_vec());
    }
    #[test]
    fn provider_rejects_unknown_contig() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, ">chr1\nACGT").unwrap();
        let p = FastaSequenceProvider::from_file(f.path()).unwrap();
        assert!(p.sequence("missing", 0, 1).is_err());
    }
}
