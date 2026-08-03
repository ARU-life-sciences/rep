use std::path::PathBuf;

use crate::{Error, ErrorKind, Result};
use csv::ReaderBuilder;

// use the CSV reader builder to parse
// blast outfmt 7, which includes headers
// we want a struct to wrap the data
// and implementations to parse
// column headers are: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

#[derive(Clone)]
pub struct BlastRecord {
    pub qseqid: String,
    pub sseqid: String,
    pub pident: f64,
    pub length: u64,
    pub mismatch: u64,
    pub gapopen: u64,
    pub qstart: u64,
    pub qend: u64,
    pub sstart: u64,
    pub send: u64,
    pub evalue: f64,
    pub bitscore: f64,
}

impl BlastRecord {
    fn from_row(row: &csv::StringRecord) -> Result<Self> {
        if row.len() != 12 {
            return Err(Error::new(ErrorKind::GenericCli(format!(
                "Expected 12 BLAST columns, found {} in record {:?}",
                row.len(),
                row
            ))));
        }

        Ok(BlastRecord {
            qseqid: row[0].to_string(),
            sseqid: row[1].to_string(),
            pident: row[2].parse()?,
            length: row[3].parse()?,
            mismatch: row[4].parse()?,
            gapopen: row[5].parse()?,
            qstart: row[6].parse()?,
            qend: row[7].parse()?,
            sstart: row[8].parse()?,
            send: row[9].parse()?,
            evalue: row[10].parse()?,
            bitscore: row[11].parse()?,
        })
    }

    pub fn from_file(path: PathBuf) -> Result<BlastTable> {
        let mut rdr = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .comment(Some(b'#'))
            .flexible(true)
            .from_path(path)?;

        let mut records = Vec::new();
        for result in rdr.records() {
            let record = Self::from_row(&result?)?;
            records.push(record);
        }

        Ok(BlastTable(records))
    }
}

#[derive(Clone)]
pub struct BlastTable(pub Vec<BlastRecord>);

impl BlastTable {
    pub fn sort_by_alignment_positions(&mut self) {
        self.0
            .sort_by(|a, b| a.sstart.partial_cmp(&b.sstart).unwrap());
    }

    pub fn filter_by_query_subject(&self, query: &str, subject: &str) -> Self {
        let filtered: Vec<BlastRecord> = self
            .0
            .iter()
            .filter(|x| x.qseqid == query && x.sseqid == subject)
            .cloned()
            .collect();
        BlastTable(filtered)
    }

    pub fn filter_by_query_name(&self, query: &str) -> Self {
        let filtered: Vec<BlastRecord> = self
            .0
            .iter()
            .filter(|x| x.qseqid == query)
            .cloned()
            .collect();
        BlastTable(filtered)
    }

    pub fn sort_by_evalue(&mut self) {
        // FIXME: sort out this unwrap
        self.0
            .sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap());
    }

    pub fn sort_by_sequence_start(&mut self) {
        self.0
            .sort_by(|a, b| a.sstart.partial_cmp(&b.sstart).unwrap());
    }

    pub fn top_n(self, n: usize) -> Self {
        BlastTable(self.0.into_iter().take(n).collect())
    }

    // Extract unique combinations of query and subject
    pub fn filter_unique_combinations(self) -> Self {
        let mut unique: Vec<BlastRecord> = Vec::new();
        for record in self.0 {
            if !unique
                .iter()
                .any(|x| x.qseqid == record.qseqid && x.sseqid == record.sseqid)
            {
                unique.push(record);
            }
        }
        BlastTable(unique)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn parses_commented_outfmt_7() {
        let directory = tempdir().unwrap();
        let path = directory.path().join("blast.tsv");
        fs::write(
            &path,
            "# BLASTN 2.15.0+\n# Fields: query id, subject id\nquery1\tchr1\t100.0\t29\t0\t0\t1\t29\t19\t47\t4.6e-14\t54.7\n# 1 hits found\n",
        )
        .unwrap();

        let table = BlastRecord::from_file(path).unwrap();
        assert_eq!(table.0.len(), 1);
        assert_eq!(table.0[0].qseqid, "query1");
        assert_eq!(table.0[0].sseqid, "chr1");
    }

    #[test]
    fn rejects_short_rows_without_panicking() {
        let directory = tempdir().unwrap();
        let path = directory.path().join("blast.tsv");
        fs::write(&path, "query1\tchr1\t100.0\n").unwrap();

        let error = match BlastRecord::from_file(path) {
            Ok(_) => panic!("short BLAST row unexpectedly parsed"),
            Err(error) => error,
        };
        assert!(error.to_string().contains("Expected 12 BLAST columns"));
    }
}
