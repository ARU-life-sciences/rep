//! A layered, versioned TE reference database.
//!
//! ```text
//! <db>/
//! ├── refdb.json                  manifest: registered source releases and last build
//! ├── sources/<name>/<version>/   immutable downloads, checksum-verified
//! ├── local/records.jsonl         our own families (append-only, with provenance)
//! └── build/                      merged libraries used by `rep classify`
//!     ├── te_nucl.fa              nucleotide consensi, `id#Class/Superfamily`
//!     ├── te_prot.fa              TE proteins, `id#Class/Superfamily`
//!     └── records.tsv             every input record, incl. duplicates and renames
//! ```
//!
//! Public sources are Repbase (EMBL releases, CC0 since 2026-06-25, mirrored by
//! Dfam) and Dfam (curated RepeatMasker library and RepeatPeps). The local layer
//! holds families found by this lab, so they are reused in later analyses.

mod embl;
mod fetch;
mod hmms;
mod labels;
mod mosaics;

pub use embl::{parse_embl, EmblEntry};
pub use fetch::{fetch, register_source, FetchOptions};
pub use fetch::{DFAM_CURRENT_URL, REPBASE_ARCHIVE_URL, REPBASE_MANIFEST_URL};
pub use hmms::{build_te_hmms, build_te_hmms_from, HmmBuildInfo, HmmBuildOptions};
pub use labels::repbase_label;
pub use mosaics::{
    find_mosaics, foreign_class, mosaics_from_hits, read_mosaics, Mosaic, MosaicInfo, MosaicOptions,
};

use std::{
    collections::{BTreeMap, HashMap},
    fmt::Write as _,
    fs,
    io::{BufRead, BufReader, Write},
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

use bio::io::fasta;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::{classify::split_label, Error, ErrorKind, Result};

pub const SCHEMA_VERSION: &str = "rep.refdb.v1";
pub const REPBASE: &str = "repbase";
pub const DFAM_CURATED: &str = "dfam-curated";
pub const DFAM_PEPTIDES: &str = "dfam-peptides";

#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum SeqKind {
    Nucleotide,
    Protein,
}

impl SeqKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Nucleotide => "nucleotide",
            Self::Protein => "protein",
        }
    }
}

impl std::str::FromStr for SeqKind {
    type Err = String;
    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        match value {
            "nucleotide" | "nucl" | "dna" => Ok(Self::Nucleotide),
            "protein" | "prot" => Ok(Self::Protein),
            _ => Err(format!("unknown sequence kind: {value}")),
        }
    }
}

/// One registered release of a public source.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct SourceRelease {
    pub name: String,
    pub version: String,
    pub url: Option<String>,
    /// Paths relative to the database root.
    pub files: Vec<String>,
    pub sha256: Vec<String>,
    pub licence: String,
    pub registered_at: u64,
}

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
#[serde(default)]
pub struct BuildInfo {
    pub built_at: u64,
    pub sources: Vec<String>,
    pub records_by_source: BTreeMap<String, usize>,
    pub nucleotide_records: usize,
    pub protein_records: usize,
    pub duplicates: usize,
    pub renamed: usize,
    pub unlabelled: usize,
    pub fingerprint: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct DbManifest {
    pub schema_version: String,
    pub sources: Vec<SourceRelease>,
    pub last_build: Option<BuildInfo>,
}

impl Default for DbManifest {
    fn default() -> Self {
        Self {
            schema_version: SCHEMA_VERSION.into(),
            sources: Vec::new(),
            last_build: None,
        }
    }
}

impl DbManifest {
    /// The most recently registered release of each source.
    pub fn current_sources(&self) -> BTreeMap<&str, &SourceRelease> {
        let mut current = BTreeMap::new();
        for release in &self.sources {
            current.insert(release.name.as_str(), release);
        }
        current
    }
}

/// A family in the local layer.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct LocalRecord {
    pub id: String,
    pub kind: SeqKind,
    pub label: String,
    /// `accepted`, `candidate` or `rejected`; rejected records are not built.
    pub status: String,
    /// Where the family came from, e.g. a genome or project name.
    pub origin: String,
    pub source_id: String,
    pub added_at: u64,
    pub notes: String,
    pub sha256: String,
    pub seq: String,
}

/// Paths of a reference database.
#[derive(Clone, Debug)]
pub struct RefDb {
    pub root: PathBuf,
}

impl RefDb {
    pub fn new(root: impl Into<PathBuf>) -> Self {
        Self { root: root.into() }
    }
    pub fn manifest_path(&self) -> PathBuf {
        self.root.join("refdb.json")
    }
    pub fn local_path(&self) -> PathBuf {
        self.root.join("local").join("records.jsonl")
    }
    pub fn build_dir(&self) -> PathBuf {
        self.root.join("build")
    }
    pub fn nucleotide_library(&self) -> PathBuf {
        self.build_dir().join("te_nucl.fa")
    }
    pub fn protein_library(&self) -> PathBuf {
        self.build_dir().join("te_prot.fa")
    }

    pub fn load_manifest(&self) -> Result<DbManifest> {
        match fs::read(self.manifest_path()) {
            Ok(bytes) => serde_json::from_slice(&bytes)
                .map_err(|e| cli_error(format!("invalid {}: {e}", self.manifest_path().display()))),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(DbManifest::default()),
            Err(e) => Err(e.into()),
        }
    }

    pub fn save_manifest(&self, manifest: &DbManifest) -> Result<()> {
        fs::create_dir_all(&self.root)?;
        let tmp = self.manifest_path().with_extension("json.tmp");
        fs::write(
            &tmp,
            serde_json::to_vec_pretty(manifest).map_err(|e| cli_error(e.to_string()))?,
        )?;
        fs::rename(tmp, self.manifest_path())?;
        Ok(())
    }

    pub fn load_local(&self) -> Result<Vec<LocalRecord>> {
        let file = match fs::File::open(self.local_path()) {
            Ok(file) => file,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(Vec::new()),
            Err(e) => return Err(e.into()),
        };
        let mut records = Vec::new();
        for (number, line) in BufReader::new(file).lines().enumerate() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            records.push(serde_json::from_str(&line).map_err(|e| {
                cli_error(format!(
                    "{} line {}: {e}",
                    self.local_path().display(),
                    number + 1
                ))
            })?);
        }
        Ok(records)
    }

    /// Append families to the local layer. Records whose ID and sequence are
    /// already present are skipped, so re-running an import is harmless; an
    /// existing ID with a different sequence is an error.
    pub fn add_local(&self, records: Vec<LocalRecord>) -> Result<AddSummary> {
        let existing = self.load_local()?;
        let mut by_id: HashMap<String, String> = existing
            .iter()
            .map(|r| (r.id.clone(), r.sha256.clone()))
            .collect();
        let mut summary = AddSummary::default();
        let mut new = Vec::new();
        for record in records {
            match by_id.get(&record.id) {
                Some(sha) if *sha == record.sha256 => summary.unchanged += 1,
                Some(_) => {
                    return Err(cli_error(format!(
                        "local record {} already exists with a different sequence; \
                         choose another --origin or retire the old record",
                        record.id
                    )))
                }
                None => {
                    by_id.insert(record.id.clone(), record.sha256.clone());
                    new.push(record);
                }
            }
        }
        if !new.is_empty() {
            fs::create_dir_all(self.local_path().parent().expect("has parent"))?;
            let mut file = fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(self.local_path())?;
            for record in &new {
                let line = serde_json::to_string(record).map_err(|e| cli_error(e.to_string()))?;
                writeln!(file, "{line}")?;
            }
        }
        summary.added = new.len();
        Ok(summary)
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AddSummary {
    pub added: usize,
    pub unchanged: usize,
}

/// Build local records from a FASTA file. Labels come from `name#label`
/// identifiers; unlabelled records become `Unknown`. A non-empty `origin` is
/// prefixed to IDs so that, e.g., `rnd-1_family-1` from two genomes stay distinct.
pub fn local_records_from_fasta(
    path: &Path,
    kind: SeqKind,
    origin: &str,
    status: &str,
    notes: &str,
) -> Result<Vec<LocalRecord>> {
    validate_status(status)?;
    let mut records = Vec::new();
    for record in fasta::Reader::from_file(path)?.records() {
        let record = record?;
        let (name, label) = split_label(record.id());
        records.push(local_record(
            &name,
            kind,
            label
                .as_deref()
                .map(normalize_label)
                .as_deref()
                .unwrap_or("Unknown"),
            origin,
            status,
            record
                .desc()
                .map_or_else(|| notes.to_owned(), |d| join_notes(notes, d)),
            record.seq(),
        ));
    }
    Ok(records)
}

/// Build local records from a `rep classify` output directory, taking families
/// whose disposition is in `dispositions` (e.g. `novel_candidate`).
pub fn local_records_from_classify(
    dir: &Path,
    dispositions: &[String],
    origin: &str,
    status: &str,
) -> Result<Vec<LocalRecord>> {
    validate_status(status)?;
    let mut sequences = HashMap::new();
    for record in fasta::Reader::from_file(dir.join("queries.fa"))?.records() {
        let record = record?;
        sequences.insert(record.id().to_owned(), record.seq().to_vec());
    }
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(dir.join("classification.tsv"))?;
    let headers = reader.headers()?.clone();
    let column = |name: &str| {
        headers
            .iter()
            .position(|h| h == name)
            .ok_or_else(|| cli_error(format!("classification.tsv has no {name} column")))
    };
    let (id_col, disp_col, label_col, conf_col, ev_col) = (
        column("family_id")?,
        column("disposition")?,
        column("label")?,
        column("confidence")?,
        column("evidence")?,
    );
    let mut records = Vec::new();
    for row in reader.records() {
        let row = row?;
        if !dispositions.iter().any(|d| d == &row[disp_col]) {
            continue;
        }
        let seq = sequences
            .get(&row[id_col])
            .ok_or_else(|| cli_error(format!("{} missing from queries.fa", &row[id_col])))?;
        let (name, _) = split_label(&row[id_col]);
        records.push(local_record(
            &name,
            SeqKind::Nucleotide,
            &row[label_col],
            origin,
            status,
            format!(
                "rep classify: {} ({} confidence); {}",
                &row[disp_col], &row[conf_col], &row[ev_col]
            ),
            seq,
        ));
    }
    Ok(records)
}

fn local_record(
    name: &str,
    kind: SeqKind,
    label: &str,
    origin: &str,
    status: &str,
    notes: String,
    seq: &[u8],
) -> LocalRecord {
    let seq = normalize_seq(seq, kind);
    let id = if origin.is_empty() {
        sanitize_id(name)
    } else {
        sanitize_id(&format!("{origin}_{name}"))
    };
    LocalRecord {
        id,
        kind,
        label: label.to_owned(),
        status: status.to_owned(),
        origin: origin.to_owned(),
        source_id: name.to_owned(),
        added_at: now(),
        notes,
        sha256: sha256_hex(seq.as_bytes()),
        seq,
    }
}

fn validate_status(status: &str) -> Result<()> {
    if matches!(status, "accepted" | "candidate" | "rejected") {
        Ok(())
    } else {
        Err(cli_error(format!(
            "status must be accepted, candidate or rejected, not {status}"
        )))
    }
}

/// One record entering a build, from any source.
#[derive(Clone, Debug)]
struct BuildRecord {
    id: String,
    kind: SeqKind,
    label: Option<String>,
    source: String,
    status: String,
    source_id: String,
    species: Option<String>,
    seq: String,
}

/// Merge the current release of each source and the local layer into the
/// build libraries. Priority for deduplication and ID ownership: accepted
/// local records, Dfam curated, Dfam peptides, Repbase, candidate local records.
pub fn build(db: &RefDb) -> Result<BuildInfo> {
    let mut manifest = db.load_manifest()?;
    let current = manifest.current_sources();
    let local = db.load_local()?;
    let mut inputs: Vec<BuildRecord> = Vec::new();
    let mut used_sources = Vec::new();

    let local_with = |status: &str| {
        local
            .iter()
            .filter(|r| r.status == status)
            .map(|r| BuildRecord {
                id: r.id.clone(),
                kind: r.kind,
                label: Some(r.label.clone()),
                source: "local".into(),
                status: r.status.clone(),
                source_id: r.source_id.clone(),
                species: Some(r.origin.clone()).filter(|o| !o.is_empty()),
                seq: r.seq.clone(),
            })
            .collect::<Vec<_>>()
    };
    inputs.extend(local_with("accepted"));
    for name in [DFAM_CURATED, DFAM_PEPTIDES, REPBASE] {
        let Some(release) = current.get(name) else {
            continue;
        };
        used_sources.push(format!("{name}@{}", release.version));
        let tag = format!("{name}@{}", release.version);
        let files = release
            .files
            .iter()
            .map(|f| db.root.join(f))
            .collect::<Vec<_>>();
        match name {
            REPBASE => inputs.extend(read_repbase(&files, &tag)?),
            DFAM_CURATED => inputs.extend(read_labelled_fasta(&files, SeqKind::Nucleotide, &tag)?),
            _ => inputs.extend(read_labelled_fasta(&files, SeqKind::Protein, &tag)?),
        }
    }
    if !local.is_empty() {
        used_sources.push(format!("local@{}", local.len()));
    }
    inputs.extend(local_with("candidate"));

    // deduplicate on (kind, sequence) and give every kept record a unique ID
    let mut seen: HashMap<(SeqKind, String), String> = HashMap::new();
    let mut ids: HashMap<String, usize> = HashMap::new();
    let mut info = BuildInfo {
        built_at: now(),
        sources: used_sources,
        ..Default::default()
    };
    let tmp = db.root.join("build.tmp");
    if tmp.exists() {
        fs::remove_dir_all(&tmp)?;
    }
    fs::create_dir_all(&tmp)?;
    let mut nucl = std::io::BufWriter::new(fs::File::create(tmp.join("te_nucl.fa"))?);
    let mut prot = std::io::BufWriter::new(fs::File::create(tmp.join("te_prot.fa"))?);
    let mut table = std::io::BufWriter::new(fs::File::create(tmp.join("records.tsv"))?);
    writeln!(
        table,
        "id\tkind\tlabel\tsource\tstatus\tsource_id\tspecies\tlength\tsha256\tduplicate_of"
    )?;
    let mut hasher_all = Sha256::new();
    for record in inputs {
        let sha = sha256_hex(record.seq.as_bytes());
        let label = record.label.clone().unwrap_or_else(|| "Unknown".into());
        let key = (record.kind, sha.clone());
        let duplicate_of = seen.get(&key).cloned();
        let mut id = record.id.clone();
        if duplicate_of.is_none() {
            let n = ids.entry(id.clone()).or_insert(0);
            *n += 1;
            if *n > 1 {
                id = format!("{}__{}", record.id, *n);
                info.renamed += 1;
            }
            seen.insert(key, id.clone());
            if record.label.is_none() {
                info.unlabelled += 1;
            }
            *info
                .records_by_source
                .entry(record.source.clone())
                .or_default() += 1;
            let header = format!(
                ">{id}#{label} source={}{}{}",
                record.source,
                if record.status == "candidate" {
                    " status=candidate"
                } else {
                    ""
                },
                record
                    .species
                    .as_deref()
                    .map(|s| format!(" [{s}]"))
                    .unwrap_or_default(),
            );
            let out = match record.kind {
                SeqKind::Nucleotide => {
                    info.nucleotide_records += 1;
                    &mut nucl
                }
                SeqKind::Protein => {
                    info.protein_records += 1;
                    &mut prot
                }
            };
            writeln!(out, "{header}")?;
            for chunk in record.seq.as_bytes().chunks(60) {
                out.write_all(chunk)?;
                out.write_all(b"\n")?;
            }
            hasher_all.update(header.as_bytes());
            hasher_all.update(record.seq.as_bytes());
        } else {
            info.duplicates += 1;
        }
        writeln!(
            table,
            "{id}\t{}\t{label}\t{}\t{}\t{}\t{}\t{}\t{sha}\t{}",
            record.kind.as_str(),
            record.source,
            record.status,
            record.source_id,
            record.species.as_deref().unwrap_or(""),
            record.seq.len(),
            duplicate_of.as_deref().unwrap_or(""),
        )?;
    }
    for mut writer in [nucl, prot, table] {
        writer.flush()?;
    }
    info.fingerprint = format!("{:x}", hasher_all.finalize());
    fs::write(
        tmp.join("build.json"),
        serde_json::to_vec_pretty(&info).map_err(|e| cli_error(e.to_string()))?,
    )?;
    if db.build_dir().exists() {
        fs::remove_dir_all(db.build_dir())?;
    }
    fs::rename(&tmp, db.build_dir())?;
    manifest.last_build = Some(info.clone());
    db.save_manifest(&manifest)?;
    Ok(info)
}

fn read_repbase(files: &[PathBuf], tag: &str) -> Result<Vec<BuildRecord>> {
    let mut records = Vec::new();
    for root in files {
        let mut paths = walkdir::WalkDir::new(root)
            .into_iter()
            .filter_map(|e| e.ok())
            .filter(|e| e.file_type().is_file())
            .map(|e| e.into_path())
            .filter(|p| p.extension().is_some_and(|x| x == "ref" || x == "embl"))
            .collect::<Vec<_>>();
        // main collection before appendix, in a stable order
        paths.sort_by_key(|p| (p.to_string_lossy().contains("appendix"), p.clone()));
        for path in paths {
            let entries = parse_embl(BufReader::new(fs::File::open(&path)?))
                .map_err(|e| cli_error(format!("{}: {e}", path.display())))?;
            for entry in entries {
                let label = repbase_label(&entry.keywords);
                for (i, protein) in entry.translations.iter().enumerate() {
                    records.push(BuildRecord {
                        id: format!("{}_{}p", entry.id, i + 1),
                        kind: SeqKind::Protein,
                        label: label.clone(),
                        source: tag.into(),
                        status: "reference".into(),
                        source_id: entry.id.clone(),
                        species: entry.species.clone(),
                        seq: normalize_seq(protein.as_bytes(), SeqKind::Protein),
                    });
                }
                records.push(BuildRecord {
                    id: entry.id.clone(),
                    kind: SeqKind::Nucleotide,
                    label,
                    source: tag.into(),
                    status: "reference".into(),
                    source_id: entry.id,
                    species: entry.species,
                    seq: normalize_seq(&entry.seq, SeqKind::Nucleotide),
                });
            }
        }
    }
    Ok(records)
}

/// FASTA whose identifiers are `name#label` (Dfam RepeatMasker libraries, RepeatPeps).
fn read_labelled_fasta(files: &[PathBuf], kind: SeqKind, tag: &str) -> Result<Vec<BuildRecord>> {
    let mut records = Vec::new();
    for path in files {
        for record in fasta::Reader::from_file(path)?.records() {
            let record = record?;
            let (name, label) = split_label(record.id());
            // Dfam libraries give the taxon as `@Clade`
            let species = record.desc().and_then(|d| {
                d.split_whitespace()
                    .find_map(|w| w.strip_prefix('@'))
                    .map(str::to_owned)
            });
            records.push(BuildRecord {
                id: name.clone(),
                kind,
                label: label.as_deref().map(normalize_label),
                source: tag.into(),
                status: "reference".into(),
                source_id: name,
                species,
                seq: normalize_seq(record.seq(), kind),
            });
        }
    }
    Ok(records)
}

/// Harmonise legacy labels with current RepeatMasker naming.
pub fn normalize_label(label: &str) -> String {
    match label {
        "LINE/Penelope" => "PLE/Penelope".into(),
        other => other.to_owned(),
    }
}

fn normalize_seq(seq: &[u8], kind: SeqKind) -> String {
    let seq = seq
        .iter()
        .filter(|b| b.is_ascii_alphabetic() || **b == b'*')
        .map(u8::to_ascii_uppercase);
    match kind {
        SeqKind::Nucleotide => seq.map(char::from).collect(),
        // terminal stops are not part of the protein
        SeqKind::Protein => seq
            .map(char::from)
            .collect::<String>()
            .trim_end_matches('*')
            .to_owned(),
    }
}

fn sanitize_id(id: &str) -> String {
    id.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.') {
                c
            } else {
                '_'
            }
        })
        .collect()
}

fn join_notes(a: &str, b: &str) -> String {
    match (a.is_empty(), b.is_empty()) {
        (true, _) => b.to_owned(),
        (_, true) => a.to_owned(),
        _ => format!("{a}; {b}"),
    }
}

pub(crate) fn sha256_hex(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

pub(crate) fn sha256_file(path: &Path) -> Result<String> {
    let mut hasher = Sha256::new();
    std::io::copy(&mut fs::File::open(path)?, &mut hasher)?;
    Ok(format!("{:x}", hasher.finalize()))
}

pub(crate) fn now() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs()
}

pub(crate) fn cli_error(message: impl Into<String>) -> Error {
    Error::new(ErrorKind::GenericCli(message.into()))
}

/// Human-readable status of a database.
pub fn status_report(db: &RefDb) -> Result<String> {
    let manifest = db.load_manifest()?;
    let mut text = String::new();
    writeln!(text, "Reference database: {}", db.root.display()).ok();
    if manifest.sources.is_empty() {
        writeln!(
            text,
            "  no public sources registered (run `rep refdb fetch`)"
        )
        .ok();
    }
    for (name, release) in manifest.current_sources() {
        writeln!(
            text,
            "  {name:<14} {} ({}){}",
            release.version,
            release.licence,
            release
                .url
                .as_deref()
                .map(|u| format!(" <{u}>"))
                .unwrap_or_default()
        )
        .ok();
    }
    let local = db.load_local()?;
    let mut by_status: BTreeMap<&str, usize> = BTreeMap::new();
    for r in &local {
        *by_status.entry(r.status.as_str()).or_default() += 1;
    }
    writeln!(
        text,
        "  local          {} records {by_status:?}",
        local.len()
    )
    .ok();
    match &manifest.last_build {
        Some(b) => writeln!(
            text,
            "  last build     {} nucleotide, {} protein, {} duplicates, {} renamed, {} unlabelled; sources {:?}",
            b.nucleotide_records, b.protein_records, b.duplicates, b.renamed, b.unlabelled, b.sources
        ),
        None => writeln!(text, "  not built yet (run `rep refdb build`)"),
    }
    .ok();
    let hmms = fs::read(db.hmm_dir().join("te_hmms.json"))
        .ok()
        .and_then(|b| serde_json::from_slice::<HmmBuildInfo>(&b).ok());
    match (hmms, &manifest.last_build) {
        (Some(h), build) => {
            let stale = build
                .as_ref()
                .is_some_and(|b| b.fingerprint != h.source_fingerprint);
            writeln!(
                text,
                "  profile HMMs   {} models for {} labels from {} proteins{}",
                h.models,
                h.labels,
                h.proteins,
                if stale {
                    " (stale: rerun `rep refdb hmms`)"
                } else {
                    ""
                }
            )
        }
        (None, _) => writeln!(text, "  profile HMMs   none (run `rep refdb hmms`)"),
    }
    .ok();
    let mosaics = fs::read(db.build_dir().join("mosaics.json"))
        .ok()
        .and_then(|b| serde_json::from_slice::<MosaicInfo>(&b).ok());
    match (mosaics, &manifest.last_build) {
        (Some(m), build) => {
            let stale = build
                .as_ref()
                .is_some_and(|b| b.fingerprint != m.source_fingerprint);
            writeln!(
                text,
                "  mosaics        {} segments in {} references{}",
                m.segments,
                m.mosaic_references,
                if stale {
                    " (stale: rerun `rep refdb mosaics`)"
                } else {
                    ""
                }
            )
        }
        (None, _) => writeln!(
            text,
            "  mosaics        not searched (run `rep refdb mosaics`)"
        ),
    }
    .ok();
    let classifier = fs::read(db.classifier())
        .ok()
        .and_then(|b| serde_json::from_slice::<serde_json::Value>(&b).ok());
    match classifier {
        Some(c) => writeln!(
            text,
            "  classifier     {} classes, trained on {} rows ({})",
            c["classes"].as_array().map_or(0, Vec::len),
            c["training"]["rows"].as_u64().unwrap_or(0),
            c["features_schema"].as_str().unwrap_or("?")
        ),
        None => writeln!(
            text,
            "  classifier     none (fit with scripts/classifier/fit.py into build/classifier.json)"
        ),
    }
    .ok();
    Ok(text)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EMBL: &str = "ID   Gypsy-1_XX   repbase; DNA; INV; 12 BP.
KW   Gypsy; LTR Retrotransposon; Transposable Element;
OS   Testus exampli
FT                   /translation=\"MKVLAWWHHK\"
SQ   Sequence 12 BP;
     acgtacgtac gt   12
//
ID   Odd-1_XX   repbase; DNA; INV; 8 BP.
KW   Nonautonomous; Transposable Element;
SQ   Sequence 8 BP;
     ttttgggg   8
//
";

    fn register(db: &RefDb, name: &str, version: &str, file: &str, contents: &str) {
        let dir = db.root.join("sources").join(name).join(version);
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join(file), contents).unwrap();
        let mut m = db.load_manifest().unwrap();
        m.sources.push(SourceRelease {
            name: name.into(),
            version: version.into(),
            url: None,
            files: vec![format!("sources/{name}/{version}/{file}")],
            sha256: vec![],
            licence: "CC0-1.0".into(),
            registered_at: 0,
        });
        db.save_manifest(&m).unwrap();
    }

    #[test]
    fn builds_merged_deduplicated_libraries_with_local_layer() {
        let dir = tempfile::tempdir().unwrap();
        let db = RefDb::new(dir.path());
        register(&db, REPBASE, "31.09", "x.ref", EMBL);
        // Dfam carries the same Gypsy-1_XX sequence and a clashing name
        register(
            &db,
            DFAM_CURATED,
            "4.0",
            "lib.fa",
            ">Gypsy-1_XX#LTR/Gypsy @Insecta\nACGTACGTACGT\n>Odd-1_XX#DNA/hAT\nAAAACCCC\n",
        );
        let fa = dir.path().join("mine.fa");
        fs::write(
            &fa,
            ">rnd-1_family-3#LTR/Unknown\nGGGGAAAATTTT\n>x#LINE/Penelope\nCCCCAAAA\n",
        )
        .unwrap();
        let recs =
            local_records_from_fasta(&fa, SeqKind::Nucleotide, "Ec32", "candidate", "").unwrap();
        assert_eq!(recs[0].id, "Ec32_rnd-1_family-3");
        assert_eq!(recs[1].label, "PLE/Penelope");
        assert_eq!(db.add_local(recs.clone()).unwrap().added, 2);
        // idempotent re-import
        assert_eq!(
            db.add_local(recs).unwrap(),
            AddSummary {
                added: 0,
                unchanged: 2
            }
        );

        let info = build(&db).unwrap();
        assert_eq!(
            info.duplicates, 1,
            "Repbase Gypsy-1_XX duplicates the Dfam copy"
        );
        assert_eq!(
            info.renamed, 1,
            "Repbase Odd-1_XX clashes with a different Dfam sequence"
        );
        assert_eq!(info.unlabelled, 1);
        assert_eq!(info.protein_records, 1);
        let nucl = fs::read_to_string(db.nucleotide_library()).unwrap();
        assert!(nucl.contains(">Gypsy-1_XX#LTR/Gypsy source=dfam-curated@4.0 [Insecta]\n"));
        assert!(nucl.contains(">Odd-1_XX__2#Unknown source=repbase@31.09\n"));
        assert!(nucl
            .contains(">Ec32_rnd-1_family-3#LTR/Unknown source=local status=candidate [Ec32]\n"));
        let prot = fs::read_to_string(db.protein_library()).unwrap();
        assert!(prot.starts_with(">Gypsy-1_XX_1p#LTR/Gypsy source=repbase@31.09"));
        let table = fs::read_to_string(db.build_dir().join("records.tsv")).unwrap();
        assert!(
            table.contains("\tGypsy-1_XX\n"),
            "duplicate recorded with its target"
        );
        assert!(status_report(&db).unwrap().contains("31.09"));
    }

    #[test]
    fn conflicting_local_ids_are_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let db = RefDb::new(dir.path());
        let one = local_record(
            "f",
            SeqKind::Nucleotide,
            "Unknown",
            "g",
            "accepted",
            String::new(),
            b"ACGT",
        );
        let two = local_record(
            "f",
            SeqKind::Nucleotide,
            "Unknown",
            "g",
            "accepted",
            String::new(),
            b"ACGA",
        );
        db.add_local(vec![one]).unwrap();
        assert!(db.add_local(vec![two]).is_err());
    }
}
