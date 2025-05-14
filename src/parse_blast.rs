use std::path::PathBuf;

use crate::Result;
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
        let r = BlastRecord {
            // TODO: these unwraps should be fine... ideally they should
            // be handled
            qseqid: row.get(0).unwrap().to_string(),
            sseqid: row.get(1).unwrap().to_string(),
            pident: row.get(2).unwrap().parse()?,
            length: row.get(3).unwrap().parse()?,
            mismatch: row.get(4).unwrap().parse()?,
            gapopen: row.get(5).unwrap().parse()?,
            qstart: row.get(6).unwrap().parse()?,
            qend: row.get(7).unwrap().parse()?,
            sstart: row.get(8).unwrap().parse()?,
            send: row.get(9).unwrap().parse()?,
            evalue: row.get(10).unwrap().parse()?,
            bitscore: row.get(11).unwrap().parse()?,
        };

        Ok(r)
    }

    pub fn from_file(path: PathBuf) -> Result<BlastTable> {
        let mut rdr = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
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
