use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
pub const SCHEMA_VERSION: &str = "rep.copy-recovery.v1";
#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct StageState {
    pub status: String,
    pub started_at: Option<String>,
    pub completed_at: Option<String>,
    pub input_fingerprint: Option<String>,
    pub relevant_config: Option<serde_json::Value>,
    pub commands: Vec<Vec<String>>,
    pub outputs: Vec<String>,
    pub warnings: Vec<String>,
    pub error: Option<String>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Manifest {
    pub schema_version: String,
    pub rep_version: String,
    pub coordinate_convention: String,
    pub configuration: serde_json::Value,
    #[serde(default)]
    pub input_fingerprint: String,
    pub warnings: Vec<String>,
    pub status: String,
    pub families: usize,
    #[serde(default)]
    pub stages: BTreeMap<String, StageState>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FamilyManifest {
    pub schema_version: String,
    pub family_id: String,
    pub status: String,
    pub raw_hits: usize,
    pub merged_loci: usize,
    pub selected: usize,
    pub warnings: Vec<String>,
    #[serde(default)]
    pub stages: BTreeMap<String, StageState>,
}
