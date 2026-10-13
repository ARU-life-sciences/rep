use std::path::PathBuf;

use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum SearchBackend {
    Blastn,
}

impl std::str::FromStr for SearchBackend {
    type Err = String;
    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        match value {
            "blastn" => Ok(Self::Blastn),
            _ => Err(format!("unsupported search backend: {value}")),
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Aligner {
    Mafft,
}

impl std::str::FromStr for Aligner {
    type Err = String;
    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        match value {
            "mafft" => Ok(Self::Mafft),
            _ => Err(format!("unsupported aligner: {value}")),
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct CopyRecoveryConfig {
    pub genome: PathBuf,
    pub queries: PathBuf,
    pub output_dir: PathBuf,
    pub flank_length: usize,
    pub hits_per_family: usize,
    pub max_fragment_gap: usize,
    pub orientation_threshold: f64,
    pub threads: usize,
    pub search_backend: SearchBackend,
    #[serde(default = "default_blast_task")]
    pub blast_task: String,
    #[serde(default = "default_evalue")]
    pub evalue: f64,
    pub aligner: Aligner,
    #[serde(default)]
    pub force: bool,
    #[serde(default)]
    pub resume: bool,
    #[serde(default)]
    pub diversity_filter: bool,
    #[serde(default)]
    pub reputils: Option<PathBuf>,
}

/// `blastn` at E ≤ 1e-5: on *Cloeon dipterum* it found twice the copies of
/// `dc-megablast` at 1e-10 (40% more full-length), with the same TSD calls, for 30% more time.
pub fn default_blast_task() -> String {
    "blastn".into()
}

pub fn default_evalue() -> f64 {
    1e-5
}

impl CopyRecoveryConfig {
    pub fn validate(&self) -> std::result::Result<(), String> {
        if self.hits_per_family == 0 {
            return Err("hits_per_family must be greater than zero".into());
        }
        if self.threads == 0 {
            return Err("threads must be greater than zero".into());
        }
        if !(0.0..=1.0).contains(&self.orientation_threshold) {
            return Err("orientation_threshold must be between 0 and 1".into());
        }
        if !matches!(
            self.blast_task.as_str(),
            "megablast" | "dc-megablast" | "blastn" | "blastn-short"
        ) {
            return Err(format!("unsupported BLAST task: {}", self.blast_task));
        }
        if self.evalue.is_nan() || self.evalue <= 0.0 {
            return Err("evalue must be positive".into());
        }
        if self.force && self.resume {
            return Err("--force and --resume are mutually exclusive".into());
        }
        Ok(())
    }
}
