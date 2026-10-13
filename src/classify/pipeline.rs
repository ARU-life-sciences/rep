use std::{
    collections::{BTreeMap, HashMap, HashSet},
    fmt::Write as _,
    fs,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    process::Command,
};

use serde::{Deserialize, Serialize};

use super::{
    copies::{read_copy_evidence, CopyEvidence, CopyParams},
    decide::{
        apply_copy_evidence, apply_host_domains, decide, is_non_te, respect_original_label,
        DecisionThresholds, Disposition, DispositionKind,
    },
    domains::{
        extract_te_models, hmmpress_command, hmmscan_command, interpret, parse_domtblout,
        parse_host_domains, DomainEvidence, DomainHit, TE_DOMAINS,
    },
    homology::{
        blastn_library_command, blastx_command, diamond_blastx_command, diamond_makedb_command,
        makeblastdb_nucleotide_command, makeblastdb_protein_command, parse_protein_hits, summarise,
        HomologyBackend, HomologySummary, ProteinHit,
    },
    library::{read_library, LibraryFamily},
    rmout::{copy_evidence_from_rm, load_genome, read_rm_out},
    structure::{analyse, StructuralEvidence, StructureParams},
    taxonomy::canonical_class,
    taxonomy::class_from_protein_name,
};
use crate::{command_runner::run_checked, CommandRunner, Error, ErrorKind, Result};

pub const SCHEMA_VERSION: &str = "rep.classify.v1";

#[derive(Clone, Debug, Serialize)]
pub struct ClassifyConfig {
    /// RepeatModeler/RepeatMasker library (`consensi.fa.classified`).
    pub library: PathBuf,
    /// TE protein FASTA, optionally gzip-compressed (e.g. RepBase proteins).
    pub proteins: Option<PathBuf>,
    /// TE nucleotide consensi with `name#label` IDs (e.g. a refdb build).
    pub nucleotides: Option<PathBuf>,
    /// Where search indexes are cached; the output directory when `None`.
    pub index_dir: Option<PathBuf>,
    /// Full Pfam-A.hmm; TE domain models are extracted from it for ORF scans.
    pub pfam: Option<PathBuf>,
    /// Pressed superfamily profile HMMs (`rep refdb hmms`), named `<model>#<label>`.
    pub te_hmms: Option<PathBuf>,
    /// Scan ORFs of families without TE protein evidence against all of
    /// `pfam` to recognise multicopy host genes.
    pub host_domains: bool,
    /// `rep recover-copies` output directory for copy-based evidence.
    pub copies_dir: Option<PathBuf>,
    /// RepeatMasker `.out` annotation of the genome with this library; used
    /// for copy evidence in preference to `copies_dir`.
    pub rm_out: Option<PathBuf>,
    /// Genome FASTA (plain or .gz) for TSD and flank analysis with `rm_out`.
    pub genome: Option<PathBuf>,
    pub copy_params: CopyParams,
    pub output_dir: PathBuf,
    /// Classify every family, not only `Unknown`s (useful for validation).
    pub all_families: bool,
    pub threads: usize,
    pub evalue: f64,
    /// BLASTN task and E-value for the nucleotide library search.
    pub nucl_task: String,
    pub nucl_evalue: f64,
    pub backend: HomologyBackend,
    pub structure: StructureParams,
    pub thresholds: DecisionThresholds,
}

impl ClassifyConfig {
    pub fn new(library: PathBuf, proteins: PathBuf, output_dir: PathBuf) -> Self {
        Self {
            library,
            proteins: Some(proteins),
            nucleotides: None,
            index_dir: None,
            pfam: None,
            te_hmms: None,
            host_domains: false,
            copies_dir: None,
            rm_out: None,
            genome: None,
            copy_params: CopyParams::default(),
            output_dir,
            all_families: false,
            threads: 1,
            evalue: 1e-5,
            nucl_task: "blastn".into(),
            nucl_evalue: 1e-5,
            backend: HomologyBackend::Blastx,
            structure: StructureParams::default(),
            thresholds: DecisionThresholds::default(),
        }
    }

    /// Use the protein and nucleotide libraries of a built reference database,
    /// caching search indexes inside it.
    pub fn from_refdb(library: PathBuf, db: &crate::refdb::RefDb, output_dir: PathBuf) -> Self {
        let mut config = Self::new(library, db.protein_library(), output_dir);
        config.nucleotides = Some(db.nucleotide_library());
        config.index_dir = Some(db.build_dir().join("index"));
        if db.te_hmms().is_file() {
            config.te_hmms = Some(db.te_hmms());
        }
        config
    }

    pub fn validate(&self) -> std::result::Result<(), String> {
        if self.threads == 0 {
            return Err("threads must be greater than zero".into());
        }
        if self.evalue.is_nan() || self.evalue <= 0.0 {
            return Err("evalue must be positive".into());
        }
        if self.proteins.is_none() && self.nucleotides.is_none() {
            return Err("a protein and/or nucleotide reference library is required".into());
        }
        let inputs = [
            Some(("library", &self.library)),
            self.proteins.as_ref().map(|p| ("proteins", p)),
            self.nucleotides.as_ref().map(|p| ("nucleotides", p)),
            self.pfam.as_ref().map(|p| ("pfam", p)),
            self.te_hmms.as_ref().map(|p| ("te_hmms", p)),
            self.rm_out.as_ref().map(|p| ("rm_out", p)),
            self.genome.as_ref().map(|p| ("genome", p)),
        ];
        for (label, path) in inputs.into_iter().flatten() {
            if !path.is_file() {
                return Err(format!("{label} is not a regular file: {}", path.display()));
            }
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct FamilyClassification {
    pub family_id: String,
    pub name: String,
    pub original_label: Option<String>,
    pub disposition: Disposition,
    /// `agree`/`disagree` with the original label when it was not Unknown.
    pub label_agreement: Option<String>,
    pub homology: HomologySummary,
    pub nucleotide_homology: HomologySummary,
    pub domains: DomainEvidence,
    pub copies: Option<CopyEvidence>,
    pub structure: StructuralEvidence,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ClassifySummary {
    pub schema_version: String,
    pub rep_version: String,
    pub library_families: usize,
    pub classified_families: usize,
    pub by_disposition: BTreeMap<String, usize>,
    pub by_label: BTreeMap<String, usize>,
    /// Number of protein hits whose subject name matched no taxonomy rule.
    pub unclassified_subject_hits: usize,
    pub label_agreement: BTreeMap<String, usize>,
    pub commands: Vec<Vec<String>>,
}

#[derive(Clone, Debug)]
pub struct ClassifyResult {
    pub output_dir: PathBuf,
    pub families: Vec<FamilyClassification>,
    pub summary: ClassifySummary,
}

pub fn run_classify(c: &ClassifyConfig) -> Result<ClassifyResult> {
    run_classify_with_runner(c, &crate::RealCommandRunner)
}

pub fn run_classify_with_runner(
    c: &ClassifyConfig,
    runner: &dyn CommandRunner,
) -> Result<ClassifyResult> {
    c.validate().map_err(cli_error)?;
    fs::create_dir_all(&c.output_dir)?;
    let library = read_library(&c.library)?;
    let selected = library
        .iter()
        .filter(|f| c.all_families || f.is_unknown())
        .collect::<Vec<_>>();
    eprintln!(
        "Classifying {} of {} library families",
        selected.len(),
        library.len()
    );

    let mut commands = Vec::new();
    let (mut hits, mut nucl_hits) = (Vec::new(), Vec::new());
    if !selected.is_empty() {
        let queries = c.output_dir.join("queries.fa");
        write_fasta(
            &queries,
            selected
                .iter()
                .map(|f| (f.id.as_str(), None, f.seq.as_slice())),
        )?;
        if let Some(proteins) = &c.proteins {
            hits = search_library(c, &queries, proteins, false, runner, &mut commands)?;
        }
        if let Some(nucleotides) = &c.nucleotides {
            nucl_hits = search_library(c, &queries, nucleotides, true, runner, &mut commands)?;
        }
    }
    let group = |hits: &'_ [ProteinHit]| {
        let mut by_family: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (i, hit) in hits.iter().enumerate() {
            by_family.entry(hit.query_id.clone()).or_default().push(i);
        }
        by_family
    };

    // ORFs of every selected family, then their TE protein domains
    let structures = selected
        .iter()
        .map(|f| analyse(&f.seq, &c.structure))
        .collect::<Vec<_>>();
    let orfs_path = c.output_dir.join("orfs.faa");
    write_orfs(&orfs_path, &selected, &structures)?;
    let lengths = selected
        .iter()
        .map(|f| (f.id.clone(), f.seq.len()))
        .collect::<HashMap<_, _>>();
    let has_orfs = structures.iter().any(|s| !s.orfs.is_empty());
    // superfamily profile HMMs: hits count as protein-level homology. They are
    // slow, and add nothing where BLASTX already found a strong hit.
    let explained = hits
        .iter()
        .filter(|h| h.bitscore >= STRONG_PROTEIN_BITS)
        .map(|h| h.query_id.as_str())
        .collect::<HashSet<_>>();
    let unexplained = selected
        .iter()
        .zip(&structures)
        .filter(|(f, s)| !s.orfs.is_empty() && !explained.contains(f.id.as_str()))
        .map(|(f, s)| (*f, s.clone()))
        .collect::<Vec<_>>();
    if let (Some(hmms), false) = (&c.te_hmms, unexplained.is_empty()) {
        let (families, structures): (Vec<_>, Vec<_>) = unexplained.into_iter().unzip();
        let profile_orfs = c.output_dir.join("te_hmms_orfs.faa");
        write_orfs(&profile_orfs, &families, &structures)?;
        let profile_hits = scan_te_hmms(c, hmms, &profile_orfs, &lengths, runner, &mut commands)?;
        fs::remove_file(&profile_orfs).ok();
        hits.extend(profile_hits);
    }
    let (hits_by_family, nucl_by_family) = (group(&hits), group(&nucl_hits));
    let domain_hits = match &c.pfam {
        Some(pfam) if has_orfs => {
            scan_domains(c, pfam, &orfs_path, &lengths, runner, &mut commands)?
        }
        _ => Vec::new(),
    };
    let mut domains_by_family: BTreeMap<&str, Vec<&DomainHit>> = BTreeMap::new();
    for hit in &domain_hits {
        domains_by_family
            .entry(&hit.family_id)
            .or_default()
            .push(hit);
    }

    // host-gene domains on ORFs of families with no TE protein evidence
    let host_by_family = match (&c.pfam, c.host_domains) {
        (Some(pfam), true) => {
            let te_protein = hits
                .iter()
                .filter(|h| h.class.as_deref().is_some_and(|l| !is_non_te(l)))
                .map(|h| h.query_id.as_str())
                .chain(domains_by_family.keys().copied())
                .collect::<HashSet<_>>();
            let (families, structures): (Vec<_>, Vec<_>) = selected
                .iter()
                .zip(&structures)
                .filter(|(f, s)| !s.orfs.is_empty() && !te_protein.contains(f.id.as_str()))
                .map(|(f, s)| (*f, s.clone()))
                .unzip();
            if families.is_empty() {
                BTreeMap::new()
            } else {
                scan_host_domains(c, pfam, &families, &structures, runner, &mut commands)?
            }
        }
        _ => BTreeMap::new(),
    };

    // copy evidence from an existing RepeatMasker annotation, if given
    let rm_hits = match &c.rm_out {
        Some(path) => {
            eprintln!("Reading RepeatMasker annotation {}", path.display());
            Some(read_rm_out(path)?)
        }
        None => None,
    };
    let genome = match (&rm_hits, &c.genome) {
        (Some(_), Some(path)) => {
            eprintln!("Loading genome {}", path.display());
            Some(load_genome(path)?)
        }
        _ => None,
    };

    let mut families = Vec::new();
    for (family, structure) in selected.iter().zip(structures) {
        let homology = summarise(
            &family_hits(&hits, &hits_by_family, &family.id),
            c.thresholds.min_bitscore,
        );
        let nucleotide_homology = summarise(
            &family_hits(&nucl_hits, &nucl_by_family, &family.id),
            c.thresholds.nucl_min_bitscore,
        );
        let domains = interpret(
            domains_by_family
                .get(family.id.as_str())
                .map_or(&[][..], Vec::as_slice),
        );
        let mut disposition = decide(
            &homology,
            &nucleotide_homology,
            &domains,
            &structure,
            &c.thresholds,
        );
        let copies = match (&rm_hits, &c.copies_dir) {
            (Some(rm), _) => Some(copy_evidence_from_rm(
                rm.get(&family.name).map_or(&[][..], Vec::as_slice),
                family.seq.len(),
                genome.as_ref(),
                &c.copy_params,
            )),
            (None, Some(dir)) => read_copy_evidence(
                &dir.join("families")
                    .join(crate::copy_recovery::family_dir_name(&family.id)),
                &c.copy_params,
            )?,
            (None, None) => None,
        };
        if let Some(evidence) = &copies {
            disposition = apply_copy_evidence(disposition, evidence, family.seq.len());
        }
        if let Some(host) = host_by_family.get(&family.id) {
            disposition = apply_host_domains(disposition, host, &structure, copies.as_ref());
        }
        let te_coding_evidence = homology.classes.iter().any(|c| !is_non_te(&c.class))
            || domains.inferred.is_some()
            || domains.novel.is_some();
        disposition =
            respect_original_label(disposition, family.label.as_deref(), te_coding_evidence);
        let label_agreement = family
            .label
            .as_deref()
            .filter(|_| !family.is_unknown())
            .map(|original| {
                if canonical_class(original) == canonical_class(&disposition.label) {
                    "agree".to_owned()
                } else {
                    "disagree".to_owned()
                }
            });
        families.push(FamilyClassification {
            family_id: family.id.clone(),
            name: family.name.clone(),
            original_label: family.label.clone(),
            disposition,
            label_agreement,
            homology,
            nucleotide_homology,
            domains,
            copies,
            structure,
        });
    }

    hits.append(&mut nucl_hits);
    write_outputs(c, &library, &families, &hits)?;
    let summary = summarise_run(&library, &families, &hits, commands);
    fs::write(
        c.output_dir.join("summary.json"),
        serde_json::to_vec_pretty(&summary).map_err(|e| cli_error(e.to_string()))?,
    )?;
    Ok(ClassifyResult {
        output_dir: c.output_dir.clone(),
        families,
        summary,
    })
}

fn family_hits<'a>(
    all: &'a [ProteinHit],
    index: &BTreeMap<String, Vec<usize>>,
    id: &str,
) -> Vec<&'a ProteinHit> {
    index
        .get(id)
        .map(|v| v.iter().map(|&i| &all[i]).collect())
        .unwrap_or_default()
}

/// Search the queries against a protein (BLASTX/DIAMOND) or nucleotide
/// (BLASTN) reference library. Indexes are reused while the library file is
/// unchanged.
fn search_library(
    c: &ClassifyConfig,
    queries: &Path,
    library: &Path,
    nucleotide: bool,
    runner: &dyn CommandRunner,
    commands: &mut Vec<Vec<String>>,
) -> Result<Vec<ProteinHit>> {
    let exe = crate::copy_recovery::executable;
    let mut run = |mut cmd: Command, what: &str| -> Result<std::process::Output> {
        commands.push(command_line(&cmd));
        run_checked(runner, &mut cmd, what)
    };
    let kind = if nucleotide {
        "nucleotides"
    } else {
        "proteins"
    };
    let library = if library.extension().is_some_and(|e| e == "gz") {
        let mut cmd = Command::new("gzip");
        cmd.arg("-dc").arg(library);
        let out = run(cmd, "gzip")?;
        let path = c.output_dir.join(format!("{kind}.fa"));
        fs::write(&path, out.stdout)?;
        path
    } else {
        library.to_path_buf()
    };
    let index_dir = c.index_dir.clone().unwrap_or_else(|| c.output_dir.clone());
    fs::create_dir_all(&index_dir)?;
    let backend = if nucleotide {
        "blastn"
    } else {
        match c.backend {
            HomologyBackend::Blastx => "blastx",
            HomologyBackend::Diamond => "diamond",
        }
    };
    let stamp = library_stamp(&library)?;
    let db = ensure_index(
        &index_dir,
        &format!("{kind}.{backend}"),
        "db",
        &stamp,
        |db| {
            let cmd = match (nucleotide, c.backend) {
                (true, _) => makeblastdb_nucleotide_command(exe("makeblastdb"), &library, db),
                (false, HomologyBackend::Blastx) => {
                    makeblastdb_protein_command(exe("makeblastdb"), &library, db)
                }
                (false, HomologyBackend::Diamond) => {
                    diamond_makedb_command(exe("diamond"), &library, db)
                }
            };
            run(cmd, "index construction").map(drop)
        },
    )?;
    let hits_path = c.output_dir.join(format!("{kind}_hits.tsv"));
    let cmd = match (nucleotide, c.backend) {
        (true, _) => blastn_library_command(
            exe("blastn"),
            &db,
            queries,
            &hits_path,
            c.threads,
            &c.nucl_task,
            c.nucl_evalue,
        ),
        (false, HomologyBackend::Blastx) => {
            blastx_command(exe("blastx"), &db, queries, &hits_path, c.threads, c.evalue)
        }
        (false, HomologyBackend::Diamond) => diamond_blastx_command(
            exe("diamond"),
            &db,
            queries,
            &hits_path,
            c.threads,
            c.evalue,
        ),
    };
    run(cmd, backend)?;
    if !hits_path.is_file() {
        // the search tools may skip the output file when nothing is found
        fs::write(&hits_path, "")?;
    }
    parse_protein_hits(&hits_path)
}

/// Return `<index_dir>/<name>/<file>`, building it with `build` unless the
/// directory's stamp matches. Builds happen in a private directory that is
/// then renamed into place, so concurrent runs sharing `index_dir` never see
/// a partial index. (BLAST databases cannot be renamed file by file because
/// they record their own name.)
fn ensure_index(
    index_dir: &Path,
    name: &str,
    file: &str,
    stamp: &str,
    build: impl FnOnce(&Path) -> Result<()>,
) -> Result<PathBuf> {
    let final_dir = index_dir.join(name);
    if fs::read_to_string(final_dir.join("stamp")).ok().as_deref() == Some(stamp) {
        return Ok(final_dir.join(file));
    }
    // PIDs repeat across hosts and containers, so add a timestamp
    let token = format!(
        "{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos()
    );
    let tmp = index_dir.join(format!(".tmp-{token}-{name}"));
    if tmp.exists() {
        fs::remove_dir_all(&tmp)?;
    }
    fs::create_dir_all(&tmp)?;
    if let Err(e) = build(&tmp.join(file)) {
        fs::remove_dir_all(&tmp).ok();
        return Err(e);
    }
    fs::write(tmp.join("stamp"), stamp)?;
    let current = |dir: &Path| fs::read_to_string(dir.join("stamp")).ok().as_deref() == Some(stamp);
    let mut published = false;
    for attempt in 0..5 {
        // renaming onto a missing directory is atomic; it fails if one exists
        if fs::rename(&tmp, &final_dir).is_ok() {
            published = true;
            break;
        }
        if current(&final_dir) {
            break; // another run published an up-to-date index first
        }
        // a stale index is in the way: move it aside (it may already be gone)
        let old = index_dir.join(format!(".old-{token}-{attempt}-{name}"));
        match fs::rename(&final_dir, &old) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(e.into()),
        }
    }
    if !published {
        fs::remove_dir_all(&tmp).ok();
        if !current(&final_dir) {
            return Err(cli_error(format!(
                "could not publish index {}",
                final_dir.display()
            )));
        }
    }
    // replaced indexes may still be open (NFS keeps them busy); remove what we can
    for entry in fs::read_dir(index_dir)?.flatten() {
        let n = entry.file_name().to_string_lossy().into_owned();
        if n.starts_with(".old-") && n.ends_with(&format!("-{name}")) {
            fs::remove_dir_all(entry.path()).ok();
        }
    }
    Ok(final_dir.join(file))
}

/// Identifies a library version for index reuse: path, size and modification time.
fn library_stamp(path: &Path) -> Result<String> {
    let meta = fs::metadata(path)?;
    let modified = meta
        .modified()?
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    Ok(format!(
        "{}\t{}\t{modified}\n",
        fs::canonicalize(path)?.display(),
        meta.len()
    ))
}

/// ORF proteins named `<family_id>|orf<n>|<strand>:<start>-<end>`, also used
/// for downstream annotation (hmmscan, InterProScan, ...).
/// BLASTX bit score above which a family counts as explained by protein
/// homology (≥98% agreement with RepeatClassifier in the v3 calibration).
const STRONG_PROTEIN_BITS: f64 = 200.0;

fn write_orfs(
    path: &Path,
    families: &[&LibraryFamily],
    structures: &[StructuralEvidence],
) -> Result<()> {
    let mut orfs = String::new();
    for (f, s) in families.iter().zip(structures) {
        for (i, orf) in s.orfs.iter().enumerate() {
            writeln!(
                orfs,
                ">{}|orf{}|{}:{}-{}\n{}",
                f.id,
                i + 1,
                orf.strand,
                orf.start,
                orf.end,
                orf.protein
            )
            .expect("writing to a String cannot fail");
        }
    }
    fs::write(path, orfs)?;
    Ok(())
}

/// Scan ORFs with the superfamily profile HMMs and report each domain hit as
/// a protein hit on the consensus (nucleotide coordinates, forward strand).
fn scan_te_hmms(
    c: &ClassifyConfig,
    hmms: &Path,
    orfs: &Path,
    lengths: &HashMap<String, usize>,
    runner: &dyn CommandRunner,
    commands: &mut Vec<Vec<String>>,
) -> Result<Vec<ProteinHit>> {
    let out = c.output_dir.join("te_hmms.domtblout");
    // hmmsearch reads each (large) model once; hmmscan would re-read the whole
    // pressed database for every ORF. E-values are scaled to the number of
    // models, as hmmscan would, so they do not depend on the library size.
    let models = count_models(hmms)?.max(1).to_string();
    let mut scan = Command::new(crate::copy_recovery::executable("hmmsearch"));
    scan.args([
        "-E",
        "1e-5",
        "--domE",
        "1e-5",
        "-Z",
        &models,
        "--domZ",
        &models,
        "--noali",
        "-o",
        "/dev/null",
        "--cpu",
    ])
    .arg(c.threads.to_string())
    .arg("--domtblout")
    .arg(&out)
    .arg(hmms)
    .arg(orfs);
    commands.push(command_line(&scan));
    run_checked(runner, &mut scan, "hmmsearch (TE profiles)")?;
    let mut hits = Vec::new();
    for line in fs::read_to_string(&out).unwrap_or_default().lines() {
        let f = line.split_whitespace().collect::<Vec<_>>();
        if line.starts_with('#') || f.len() < 22 {
            continue;
        }
        // hmmsearch: target (ORF) first, then query (model)
        let (orf, model, model_length) = (f[0], f[3], f[5]);
        let mut parts = orf.rsplitn(3, '|');
        let (Some(coords), Some(_), Some(family)) = (parts.next(), parts.next(), parts.next())
        else {
            continue;
        };
        let (strand, range) = coords.split_at(1);
        let Some((start, end)) = range
            .trim_start_matches(':')
            .split_once('-')
            .and_then(|(s, e)| Some((s.parse::<u64>().ok()?, e.parse::<u64>().ok()?)))
        else {
            continue;
        };
        let (from, to): (u64, u64) = (f[17].parse()?, f[18].parse()?);
        // amino-acid positions on the ORF -> 1-based nucleotide positions
        let (q_start, q_end) = if strand == "+" {
            (start + 3 * (from - 1) + 1, start + 3 * to)
        } else {
            (end - 3 * to + 1, end - 3 * (from - 1))
        };
        hits.push(ProteinHit {
            query_id: family.to_owned(),
            subject_id: model.to_owned(),
            class: class_from_protein_name(model),
            query_length: lengths.get(family).copied().unwrap_or(0) as u64,
            subject_length: model_length.parse()?,
            percent_identity: 0.0,
            query_start: q_start,
            query_end: q_end,
            strand: strand.chars().next().unwrap_or('+'),
            subject_start: f[15].parse()?,
            subject_end: f[16].parse()?,
            evalue: f[12].parse()?,
            bitscore: f[13].parse()?,
        });
    }
    Ok(hits)
}

/// Number of models in an HMM file: the `te_hmms.tsv` table written beside
/// it by `rep refdb hmms`, else the `NAME` lines of the file itself.
fn count_models(hmms: &Path) -> Result<usize> {
    let table = hmms.with_extension("tsv");
    if let Ok(text) = fs::read_to_string(&table) {
        return Ok(text.lines().skip(1).filter(|l| !l.is_empty()).count());
    }
    let mut n = 0;
    for line in BufReader::new(fs::File::open(hmms)?).lines() {
        n += usize::from(line?.starts_with("NAME "));
    }
    Ok(n)
}

/// Scan the ORFs of `families` against the full Pfam (gathering thresholds)
/// and return their host-gene domains by family.
fn scan_host_domains(
    c: &ClassifyConfig,
    pfam: &Path,
    families: &[&LibraryFamily],
    structures: &[StructuralEvidence],
    runner: &dyn CommandRunner,
    commands: &mut Vec<Vec<String>>,
) -> Result<BTreeMap<String, Vec<(String, f64)>>> {
    let orfs = c.output_dir.join("host_orfs.faa");
    write_orfs(&orfs, families, structures)?;
    let out = c.output_dir.join("host_domains.domtblout");
    let mut scan = Command::new(crate::copy_recovery::executable("hmmsearch"));
    scan.args(["--cut_ga", "--noali", "-o", "/dev/null", "--cpu"])
        .arg(c.threads.to_string())
        .arg("--domtblout")
        .arg(&out)
        .arg(pfam)
        .arg(&orfs);
    commands.push(command_line(&scan));
    run_checked(runner, &mut scan, "hmmsearch (host domains)")?;
    fs::remove_file(&orfs).ok();
    parse_host_domains(&out)
}

/// Scan ORFs with the TE domain subset of Pfam, extracting and pressing the
/// subset once per Pfam file.
fn scan_domains(
    c: &ClassifyConfig,
    pfam: &Path,
    orfs: &Path,
    lengths: &HashMap<String, usize>,
    runner: &dyn CommandRunner,
    commands: &mut Vec<Vec<String>>,
) -> Result<Vec<DomainHit>> {
    let exe = crate::copy_recovery::executable;
    let index_dir = c.index_dir.clone().unwrap_or_else(|| c.output_dir.clone());
    fs::create_dir_all(&index_dir)?;
    let stamp = format!(
        "{}{}",
        library_stamp(pfam)?,
        TE_DOMAINS
            .iter()
            .map(|(n, _)| *n)
            .collect::<Vec<_>>()
            .join(",")
    );
    let hmm = ensure_index(&index_dir, "te_domains", "te_domains.hmm", &stamp, |hmm| {
        let missing = extract_te_models(pfam, hmm)?;
        if !missing.is_empty() {
            eprintln!(
                "TE domain models absent from this Pfam: {}",
                missing.join(", ")
            );
        }
        let mut press = hmmpress_command(exe("hmmpress"), hmm);
        commands.push(command_line(&press));
        run_checked(runner, &mut press, "hmmpress").map(drop)
    })?;
    let out = c.output_dir.join("domains.domtblout");
    let mut scan = hmmscan_command(exe("hmmscan"), &hmm, orfs, &out, c.threads);
    commands.push(command_line(&scan));
    run_checked(runner, &mut scan, "hmmscan")?;
    if !out.is_file() {
        fs::write(&out, "")?;
    }
    parse_domtblout(&out, lengths)
}

fn command_line(cmd: &Command) -> Vec<String> {
    std::iter::once(cmd.get_program())
        .chain(cmd.get_args())
        .map(|s| s.to_string_lossy().into_owned())
        .collect()
}

fn write_outputs(
    c: &ClassifyConfig,
    library: &[LibraryFamily],
    families: &[FamilyClassification],
    hits: &[ProteinHit],
) -> Result<()> {
    let out = &c.output_dir;

    let mut tsv = String::from(
        "family_id\toriginal_label\tlength\tdisposition\tlabel\tconfidence\tlabel_agreement\t\
         best_class\tbest_subject\tbest_bitscore\tclass_coverage\tchimeric_classes\t\
         nucl_best_class\tnucl_best_subject\tnucl_bitscore\tnucl_coverage\tnucl_identity\tnucl_chimeric_classes\t\
         domain_architecture\tdomain_class\tdomain_novelty\t\
         copy_loci\tcopy_full_length\tcopy_identity\tgenome_bp\ttsd_length\ttsd_motif\ttsd_support\tflank_similarity\t\
         ltr_length\tltr_identity\ttir_length\ttir_identity\ttail\ttandem_period\t\
         tandem_identity\tlongest_orf_aa\tn_orfs\tgc\tevidence\twarnings\n",
    );
    for f in families {
        let best = f.homology.best();
        let nbest = f.nucleotide_homology.best();
        let s = &f.structure;
        let opt = |v: Option<String>| v.unwrap_or_default();
        writeln!(
            tsv,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}",
            f.family_id,
            f.original_label.as_deref().unwrap_or(""),
            s.length,
            f.disposition.kind.as_str(),
            f.disposition.label,
            f.disposition.confidence.as_str(),
            f.label_agreement.as_deref().unwrap_or(""),
            opt(best.map(|b| b.class.clone())),
            opt(best.map(|b| b.best_subject.clone())),
            opt(best.map(|b| format!("{:.1}", b.best_bitscore))),
            opt(best.map(|b| format!("{:.3}", b.query_coverage))),
            f.homology.chimeric_classes.join(","),
            opt(nbest.map(|b| b.class.clone())),
            opt(nbest.map(|b| b.best_subject.clone())),
            opt(nbest.map(|b| format!("{:.1}", b.best_bitscore))),
            opt(nbest.map(|b| format!("{:.3}", b.query_coverage))),
            opt(nbest.map(|b| format!("{:.1}", b.best_identity))),
            f.nucleotide_homology.chimeric_classes.join(","),
            f.domains.architecture,
            f.domains.inferred.as_deref().unwrap_or(""),
            f.domains.novel.as_deref().unwrap_or(""),
            opt(f.copies.as_ref().map(|c| c.loci.to_string())),
            opt(f.copies.as_ref().map(|c| c.full_length.to_string())),
            opt(f.copies.as_ref().map(|c| format!("{:.1}", c.mean_identity))),
            opt(f.copies.as_ref().map(|c| c.genome_bp.to_string())),
            opt(f.copies.as_ref().and_then(|c| c.tsd.as_ref()).map(|t| t.length.to_string())),
            opt(f.copies.as_ref().and_then(|c| c.tsd.as_ref()).map(|t| t.motif.clone())),
            opt(f.copies.as_ref().and_then(|c| c.tsd.as_ref()).map(|t| format!("{:.2}", t.support))),
            opt(f.copies.as_ref().and_then(|c| c.flank_similarity).map(|v| format!("{v:.2}"))),
            opt(s.ltr.as_ref().map(|r| r.length.to_string())),
            opt(s.ltr.as_ref().map(|r| format!("{:.3}", r.identity))),
            opt(s.tir.as_ref().map(|r| r.length.to_string())),
            opt(s.tir.as_ref().map(|r| format!("{:.3}", r.identity))),
            s.tail.map_or("", |t| t.as_str()),
            opt(s.tandem.as_ref().map(|t| t.period.to_string())),
            opt(s.tandem.as_ref().map(|t| format!("{:.3}", t.identity))),
            s.longest_orf_aa(),
            s.orfs.len(),
            s.composition.gc,
            f.disposition.evidence.join(";"),
            f.disposition.warnings.join(";"),
        )
        .expect("writing to a String cannot fail");
    }
    fs::write(out.join("classification.tsv"), tsv)?;
    fs::write(
        out.join("families.json"),
        serde_json::to_vec_pretty(families).map_err(|e| cli_error(e.to_string()))?,
    )?;

    // Library with Unknowns relabelled; already-classified families are kept as they are.
    let relabel = families
        .iter()
        .filter(|f| f.label_agreement.is_none() && f.disposition.label != "Unknown")
        .map(|f| (f.family_id.as_str(), f.disposition.label.as_str()))
        .collect::<BTreeMap<_, _>>();
    let renamed = library
        .iter()
        .map(|f| {
            let id = relabel
                .get(f.id.as_str())
                .map_or_else(|| f.id.clone(), |label| format!("{}#{label}", f.name));
            (id, f.description.as_deref(), f.seq.as_slice())
        })
        .collect::<Vec<_>>();
    write_fasta(
        &out.join("library.reclassified.fa"),
        renamed.iter().map(|(id, d, s)| (id.as_str(), *d, *s)),
    )?;

    let by_id = library
        .iter()
        .map(|f| (f.id.as_str(), f))
        .collect::<BTreeMap<_, _>>();
    write_fasta(
        &out.join("novel_candidates.fa"),
        families
            .iter()
            .filter(|f| f.disposition.kind == DispositionKind::NovelCandidate)
            .map(|f| {
                (
                    f.family_id.as_str(),
                    Some(f.disposition.label.as_str()),
                    by_id[f.family_id.as_str()].seq.as_slice(),
                )
            }),
    )?;

    // Protein subjects that the taxonomy could not place, to extend the rules.
    let mut unplaced: BTreeMap<&str, (usize, f64)> = BTreeMap::new();
    for hit in hits.iter().filter(|h| h.class.is_none()) {
        let e = unplaced.entry(&hit.subject_id).or_default();
        e.0 += 1;
        e.1 = e.1.max(hit.bitscore);
    }
    let mut text = String::from("subject_id\tn_hits\tbest_bitscore\n");
    for (subject, (n, bits)) in unplaced {
        writeln!(text, "{subject}\t{n}\t{bits}").expect("writing to a String cannot fail");
    }
    fs::write(out.join("unclassified_subjects.tsv"), text)?;
    Ok(())
}

fn summarise_run(
    library: &[LibraryFamily],
    families: &[FamilyClassification],
    hits: &[ProteinHit],
    commands: Vec<Vec<String>>,
) -> ClassifySummary {
    let mut by_disposition = BTreeMap::new();
    let mut by_label = BTreeMap::new();
    let mut label_agreement = BTreeMap::new();
    for f in families {
        *by_disposition
            .entry(f.disposition.kind.as_str().to_owned())
            .or_insert(0) += 1;
        *by_label.entry(f.disposition.label.clone()).or_insert(0) += 1;
        if let Some(a) = &f.label_agreement {
            *label_agreement.entry(a.clone()).or_insert(0) += 1;
        }
    }
    ClassifySummary {
        schema_version: SCHEMA_VERSION.into(),
        rep_version: env!("CARGO_PKG_VERSION").into(),
        library_families: library.len(),
        classified_families: families.len(),
        by_disposition,
        by_label,
        unclassified_subject_hits: hits.iter().filter(|h| h.class.is_none()).count(),
        label_agreement,
        commands,
    }
}

fn write_fasta<'a>(
    path: &Path,
    records: impl Iterator<Item = (&'a str, Option<&'a str>, &'a [u8])>,
) -> Result<()> {
    let mut text = String::new();
    for (id, desc, seq) in records {
        match desc {
            Some(d) => writeln!(text, ">{id} {d}"),
            None => writeln!(text, ">{id}"),
        }
        .expect("writing to a String cannot fail");
        for chunk in seq.chunks(60) {
            text.push_str(&String::from_utf8_lossy(chunk));
            text.push('\n');
        }
    }
    fs::write(path, text)?;
    Ok(())
}

fn cli_error(message: impl Into<String>) -> Error {
    Error::new(ErrorKind::GenericCli(message.into()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::process::ExitStatusExt;
    use std::process::{ExitStatus, Output};

    /// Pretends to be gzip/makeblastdb/blastx: writes canned BLASTX hits to `-out`.
    struct MockRunner {
        hits: String,
        nucl_hits: String,
    }

    impl CommandRunner for MockRunner {
        fn run(&self, cmd: &mut Command) -> Result<Output> {
            let args = cmd
                .get_args()
                .map(|a| a.to_string_lossy().into_owned())
                .collect::<Vec<_>>();
            if let Some(i) = args.iter().position(|a| a == "-out") {
                if args.iter().any(|a| a == "-query") {
                    let blastn = cmd.get_program().to_string_lossy().ends_with("blastn");
                    let hits = if blastn { &self.nucl_hits } else { &self.hits };
                    fs::write(&args[i + 1], hits)?;
                }
            }
            Ok(Output {
                status: ExitStatus::from_raw(0),
                stdout: Vec::new(),
                stderr: Vec::new(),
            })
        }
    }

    fn random(n: usize, seed: u64) -> String {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                ['A', 'C', 'G', 'T'][(x % 4) as usize]
            })
            .collect()
    }

    #[test]
    fn classifies_unknowns_and_relabels_library() {
        let dir = tempfile::tempdir().unwrap();
        let ltr = random(200, 11);
        let ltr_element = format!("{ltr}{}{ltr}", random(3000, 4));
        let library = format!(
            ">fam1#Unknown\n{}\n>fam2#Unknown\n{}\n>fam3#LINE/L1\n{}\n>fam4#Unknown\n{}\n>fam5#Unknown\n{}\n",
            random(3000, 1),
            "TTAGGG".repeat(60),
            random(2000, 3),
            ltr_element,
            random(1500, 5),
        );
        let lib = dir.path().join("consensi.fa.classified");
        fs::write(&lib, library).unwrap();
        let proteins = dir.path().join("proteins.fa");
        fs::write(&proteins, ">Gypsy-1_xx\nMKV\n").unwrap();
        let hits = "fam1#Unknown\tGypsy-1_xx\t3000\t600\t40.0\t300\t1\t1200\t1\t400\t1e-40\t350\n\
                    fam5#Unknown\tItch-1_xx\t1500\t600\t30.0\t100\t1\t300\t1\t100\t1e-10\t70\n";
        let runner = MockRunner {
            hits: hits.into(),
            nucl_hits: String::new(),
        };
        let mut config = ClassifyConfig::new(lib, proteins, dir.path().join("out"));
        config.threads = 2;
        let result = run_classify_with_runner(&config, &runner).unwrap();

        let kinds = result
            .families
            .iter()
            .map(|f| {
                (
                    f.name.as_str(),
                    f.disposition.kind,
                    f.disposition.label.as_str(),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            kinds,
            vec![
                ("fam1", DispositionKind::KnownTe, "LTR/Gypsy"),
                ("fam2", DispositionKind::SimpleRepeat, "Simple_repeat"),
                ("fam4", DispositionKind::NovelCandidate, "LTR/Unknown"),
                ("fam5", DispositionKind::Unresolved, "Unknown"),
            ]
        );
        let out = dir.path().join("out");
        let relabelled = fs::read_to_string(out.join("library.reclassified.fa")).unwrap();
        assert!(relabelled.contains(">fam1#LTR/Gypsy\n"));
        assert!(relabelled.contains(">fam3#LINE/L1\n"));
        assert!(relabelled.contains(">fam5#Unknown\n"));
        let novel = fs::read_to_string(out.join("novel_candidates.fa")).unwrap();
        assert!(novel.starts_with(">fam4#Unknown LTR/Unknown\n"));
        let unplaced = fs::read_to_string(out.join("unclassified_subjects.tsv")).unwrap();
        assert!(unplaced.contains("Itch-1_xx\t1\t70"));
        assert_eq!(
            fs::read_to_string(out.join("classification.tsv"))
                .unwrap()
                .lines()
                .count(),
            5
        );
        assert_eq!(result.summary.by_disposition["known_te"], 1);
        assert_eq!(result.summary.commands.len(), 2);
    }

    #[test]
    fn all_families_mode_reports_label_agreement() {
        let dir = tempfile::tempdir().unwrap();
        let lib = dir.path().join("lib.fa");
        fs::write(&lib, format!(">famA#LTR/Copia\n{}\n", random(3000, 9))).unwrap();
        let proteins = dir.path().join("p.fa");
        fs::write(&proteins, ">Copia-1_xx\nMKV\n").unwrap();
        let runner = MockRunner {
            hits: "famA#LTR/Copia\tCopia-1_xx\t3000\t600\t40.0\t300\t1\t1200\t1\t400\t1e-40\t350\n"
                .into(),
            nucl_hits: String::new(),
        };
        let mut config = ClassifyConfig::new(lib, proteins, dir.path().join("out"));
        config.all_families = true;
        let result = run_classify_with_runner(&config, &runner).unwrap();
        assert_eq!(result.families[0].label_agreement.as_deref(), Some("agree"));
        assert_eq!(result.summary.label_agreement["agree"], 1);
    }

    #[test]
    fn nucleotide_library_evidence_and_index_reuse() {
        let dir = tempfile::tempdir().unwrap();
        let lib = dir.path().join("lib.fa");
        fs::write(
            &lib,
            format!(
                ">sine1#Unknown\n{}\n>rrna#Unknown\n{}\n",
                random(300, 21),
                random(1800, 22)
            ),
        )
        .unwrap();
        let nucl = dir.path().join("te_nucl.fa");
        fs::write(&nucl, ">B1_Mm#SINE/7SL\nACGT\n>5S#rRNA\nACGT\n").unwrap();
        let runner = MockRunner {
            hits: String::new(),
            nucl_hits:
                "sine1#Unknown\tB1_Mm#SINE/7SL\t300\t140\t88.0\t280\t1\t280\t1\t140\t1e-60\t320\n\
                        rrna#Unknown\t5S#rRNA\t1800\t1800\t99.0\t1700\t1\t1700\t1\t1700\t0\t3000\n"
                    .into(),
        };
        let mut config = ClassifyConfig::new(lib, PathBuf::new(), dir.path().join("out"));
        config.proteins = None;
        config.nucleotides = Some(nucl);
        config.index_dir = Some(dir.path().join("index"));
        let result = run_classify_with_runner(&config, &runner).unwrap();
        let got = result
            .families
            .iter()
            .map(|f| (f.disposition.kind, f.disposition.label.as_str()))
            .collect::<Vec<_>>();
        assert_eq!(
            got,
            vec![
                (DispositionKind::KnownTe, "SINE/7SL"),
                (DispositionKind::HostOrArtefact, "rRNA"),
            ]
        );
        let builds = |r: &ClassifyResult| {
            r.summary
                .commands
                .iter()
                .filter(|c| c[0].ends_with("makeblastdb"))
                .count()
        };
        assert_eq!(builds(&result), 1);
        // the index is reused while the library is unchanged
        let again = run_classify_with_runner(&config, &runner).unwrap();
        assert_eq!(builds(&again), 0);
    }

    /// Fakes hmmsearch: one hit of a labelled profile on a forward ORF.
    struct HmmsearchRunner;

    impl CommandRunner for HmmsearchRunner {
        fn run(&self, cmd: &mut Command) -> Result<Output> {
            let args = cmd
                .get_args()
                .map(|a| a.to_string_lossy().into_owned())
                .collect::<Vec<_>>();
            let z = args.iter().position(|a| a == "-Z").unwrap();
            assert_eq!(args[z + 1], "3");
            let out = args.iter().position(|a| a == "--domtblout").unwrap();
            fs::write(
                &args[out + 1],
                "# target query\n\
                 fam#Unknown|orf1|+:100-1000 - 300 hmm_7#DNA/Crypton - 450 1e-30 100.0 0.1 1 1 \
                 1e-30 1e-30 99.0 0.1 5 205 11 210 11 210 0.9 -\n",
            )?;
            Ok(Output {
                status: ExitStatus::from_raw(0),
                stdout: Vec::new(),
                stderr: Vec::new(),
            })
        }
    }

    #[test]
    fn te_profile_hits_map_to_family_coordinates() {
        let dir = tempfile::tempdir().unwrap();
        let hmms = dir.path().join("te_hmms.hmm");
        fs::write(&hmms, "").unwrap();
        fs::write(
            dir.path().join("te_hmms.tsv"),
            "model\tlabel\tmembers\na\tx\t5\nb\tx\t5\nc\tx\t5\n",
        )
        .unwrap();
        let config = ClassifyConfig::new(PathBuf::new(), PathBuf::new(), dir.path().into());
        let lengths = HashMap::from([("fam#Unknown".to_owned(), 1200)]);
        let mut commands = Vec::new();
        let hits = scan_te_hmms(
            &config,
            &hmms,
            Path::new("orfs.faa"),
            &lengths,
            &HmmsearchRunner,
            &mut commands,
        )
        .unwrap();
        assert!(commands[0][0].ends_with("hmmsearch"));
        assert_eq!(hits.len(), 1);
        let h = &hits[0];
        assert_eq!(h.query_id, "fam#Unknown");
        assert_eq!(h.class.as_deref(), Some("DNA/Crypton"));
        assert_eq!(
            (h.query_start, h.query_end, h.query_length),
            (131, 730, 1200)
        );
        assert_eq!(
            (h.subject_start, h.subject_end, h.subject_length),
            (5, 205, 450)
        );
        assert_eq!(h.bitscore, 99.0);
    }
}
