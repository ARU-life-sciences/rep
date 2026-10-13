//! Classification and triage of repeat families, focused on RepeatModeler `Unknown`s.
//!
//! Every selected family receives a disposition backed by explicit evidence:
//! protein homology to classified TE proteins, tandem/low-complexity structure,
//! or TE hallmarks (terminal repeats, tails, coding capacity) without known
//! homology, which makes it a candidate novel element. Families that nothing
//! explains are reported as `unresolved` together with the missing evidence.

mod copies;
mod decide;
mod domains;
mod features;
mod homology;
mod library;
mod model;
mod pipeline;
mod rmout;
mod structure;
mod taxonomy;

pub use copies::{expected_tsd, read_copy_evidence, CopyEvidence, CopyParams, Tsd};
pub use decide::{
    apply_copy_evidence, apply_host_domains, decide, is_non_te, Confidence, DecisionThresholds,
    Disposition, DispositionKind, HOST_GENE_LABEL,
};
pub use domains::{
    extract_te_models, hmmscan_command, interpret, is_te_associated, parse_domtblout,
    parse_host_domains, role_of, DomainEvidence, DomainHit, Role, TE_DOMAINS,
};
pub use features::{FEATURES_SCHEMA, FEATURE_CLASSES};
pub use homology::{
    blastn_library_command, blastx_command, diamond_blastx_command, diamond_makedb_command,
    makeblastdb_nucleotide_command, makeblastdb_protein_command, parse_protein_hits, summarise,
    ClassSupport, HomologyBackend, HomologySummary, ProteinHit,
};
pub use library::{read_library, split_label, LibraryFamily};
pub use model::{Classifier, ModelCall, CLASSIFIER_SCHEMA};
pub use pipeline::{
    run_classify, run_classify_with_runner, ClassifyConfig, ClassifyResult, ClassifySummary,
    FamilyClassification, SCHEMA_VERSION,
};
pub use rmout::{copy_evidence_from_rm, load_genome, read_rm_out, RmHit};
pub use structure::{
    analyse, orfs, Orf, StructuralEvidence, StructureParams, Tail, Tandem, TerminalRepeat,
};
pub use taxonomy::{canonical_class, class_from_protein_name, top_level};
