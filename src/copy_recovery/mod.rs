//! Reusable candidate-consensus copy recovery and family-alignment pipeline.
//!
//! Search-hit coordinates are retained as BLAST's 1-based inclusive coordinates.
//! Locus and BED coordinates are normalized to 0-based half-open intervals.

mod align;
mod config;
mod extract;
mod hits;
mod loci;
mod manifest;
mod pipeline;
mod search;

pub use align::mafft_command;
pub use config::{Aligner, CopyRecoveryConfig, SearchBackend};
pub use extract::{clip_interval, reverse_complement, FastaSequenceProvider, SequenceProvider};
pub use hits::{parse_blast_file, Hit, Strand};
pub use loci::{collapse_overlapping_loci, merge_hits, rank_loci, Locus};
pub use manifest::{FamilyManifest, Manifest, StageState, SCHEMA_VERSION};
pub use pipeline::{
    run_copy_recovery, run_copy_recovery_with_runner, CopyRecoveryResult, FamilyRecoveryResult,
};
pub use search::{blastn_command, makeblastdb_command};
