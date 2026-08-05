use serde::{Deserialize, Serialize};
pub const SCHEMA_VERSION: &str = "rep.copy-recovery.v1";
#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct StageState {
    pub status: String,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Manifest {
    pub schema_version: String,
    pub rep_version: String,
    pub coordinate_convention: String,
    pub configuration: serde_json::Value,
    pub warnings: Vec<String>,
    pub status: String,
    pub families: usize,
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
}
