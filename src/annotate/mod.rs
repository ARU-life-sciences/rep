//! `rep annotate`: transposable-element fragments in any sequence (roadmap
//! M6, a CENSOR replacement).
//!
//! The input (contigs, BACs, transcripts, regions of an assembly) is cut
//! into overlapping windows and searched against the reference database and
//! any `--library` files. Bases are then claimed in four rounds, each only
//! on what the earlier rounds left:
//! 1. **Nucleotide** hits (BLASTN `-task blastn`), best bit score first.
//! 2. **Protein** hits (BLASTX against TE proteins) on unclaimed stretches,
//!    for elements too old for their DNA to match.
//! 3. **Tandem** runs (microsatellites, satellites with units up to 200 bp)
//!    and **low complexity**, from the sequence itself.
//! 4. **Structure**: terminal-repeat pairs (LTR- or TIR-like) within long
//!    unannotated stretches, reported as candidate novel elements. They are
//!    never given a class and are not masked.
//!
//! A hit that loses some of its bases keeps the rest, with its reference
//! coordinates interpolated. Fragments of the same reference on the same
//! strand, in reference order and close once insertions are discounted, are
//! chained into one element.
//!
//! Whole genomes belong to RepeatMasker in the pipeline; for a single
//! consensus, `rep classify` makes the class call.

mod resolve;
mod scan;

use std::{
    collections::{BTreeMap, HashMap},
    fmt::Write as _,
    fs,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    process::Command,
};

use bio::io::fasta;
use serde::Serialize;

use crate::{
    classify::{
        canonical_class, ensure_index, is_non_te, library_stamp, makeblastdb_nucleotide_command,
        makeblastdb_protein_command, parse_protein_hits, ProteinHit, OUTFMT_FIELDS,
    },
    command_runner::run_checked,
    copy_recovery::executable,
    refdb::RefDb,
    CommandRunner, Error, ErrorKind, Result,
};
use resolve::{chain, claim, Occupancy};
use scan::{low_complexity_runs, tandem_runs, terminal_pairs, SelfHit};

#[derive(Clone, Debug, Serialize)]
pub struct AnnotateConfig {
    pub input: PathBuf,
    /// Built `rep refdb` database: nucleotide and protein references.
    pub refdb: Option<PathBuf>,
    /// Further nucleotide libraries with `name#Class/Superfamily` headers,
    /// e.g. a classified `rep` library or a `rep relatives` file.
    pub libraries: Vec<PathBuf>,
    /// TE protein library used instead of the database's.
    pub protein_library: Option<PathBuf>,
    pub output_dir: PathBuf,
    pub threads: usize,
    pub nucl_evalue: f64,
    pub prot_evalue: f64,
    /// Weaker protein hits are mostly host genes with transposase-like
    /// domains. On D. melanogaster (dm6, 2026-10-15, Drosophila references
    /// held out) 80 bits kept 87% of the protein-found bases in the
    /// RepeatMasker track of a pericentric region and cut those outside it
    /// in a euchromatic region from 7.8 to 4.6 kb; 100 bits cut that to
    /// 2.5 kb but lost 7 points of sensitivity in the pericentric region.
    pub min_protein_bits: f64,
    /// Search window and overlap for long sequences.
    pub window: u64,
    pub window_overlap: u64,
    /// Shortest fragment kept after overlap resolution.
    pub min_fragment: u64,
    /// Largest gap (insertions discounted) bridged when chaining fragments.
    pub chain_gap: u64,
    pub proteins: bool,
    /// Shortest unclaimed stretch searched against proteins.
    pub protein_min_gap: u64,
    pub novel: bool,
    /// Shortest unannotated stretch searched for terminal-repeat pairs.
    pub novel_min_gap: u64,
    pub hard_mask: bool,
}

impl AnnotateConfig {
    pub fn new(input: PathBuf, output_dir: PathBuf) -> Self {
        Self {
            input,
            refdb: None,
            libraries: Vec::new(),
            protein_library: None,
            output_dir,
            threads: 1,
            nucl_evalue: 1e-5,
            prot_evalue: 1e-5,
            min_protein_bits: 80.0,
            window: 20_000,
            window_overlap: 2_000,
            min_fragment: 30,
            chain_gap: 500,
            proteins: true,
            protein_min_gap: 150,
            novel: true,
            novel_min_gap: 400,
            hard_mask: false,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Evidence {
    #[default]
    Nucleotide,
    Protein,
    Tandem,
    LowComplexity,
    Structure,
}

impl Evidence {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Nucleotide => "nucleotide",
            Self::Protein => "protein",
            Self::Tandem => "tandem",
            Self::LowComplexity => "low_complexity",
            Self::Structure => "structure",
        }
    }
}

/// One annotated stretch of an input sequence.
#[derive(Clone, Debug, Default, Serialize)]
pub struct Fragment {
    /// Index of the input sequence.
    pub seq: usize,
    /// 0-based half-open.
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub reference: String,
    /// `refdb`, a library file name, or `scan` for sequence-intrinsic calls.
    pub source: String,
    /// 1-based inclusive; amino acids for protein references; 0 if none.
    pub ref_start: u64,
    pub ref_end: u64,
    pub ref_length: u64,
    pub class: String,
    pub identity: f64,
    pub bitscore: f64,
    pub evalue: f64,
    pub evidence: Evidence,
    pub element: usize,
    pub note: String,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct SequenceSummary {
    pub name: String,
    pub length: u64,
    pub masked_bp: u64,
    pub by_class: BTreeMap<String, u64>,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct AnnotateSummary {
    pub sequences: usize,
    pub total_bp: u64,
    pub masked_bp: u64,
    pub fragments: usize,
    pub elements: usize,
    /// Annotated bases per top-level class (candidates excluded).
    pub by_class: BTreeMap<String, u64>,
    pub by_evidence: BTreeMap<String, u64>,
    pub candidate_novel: usize,
    pub candidate_novel_bp: u64,
    pub per_sequence: Vec<SequenceSummary>,
    /// Set when the input looks like one element rather than a region.
    pub hint: Option<String>,
    pub commands: Vec<Vec<String>>,
}

pub const TABLE: &str = "annotate.tsv";
pub const GFF: &str = "annotate.gff3";
pub const BED: &str = "annotate.bed";
pub const MASKED: &str = "masked.fa";
pub const SUMMARY_JSON: &str = "summary.json";
pub const SUMMARY_TSV: &str = "summary.tsv";

struct Input {
    names: Vec<String>,
    seqs: Vec<Vec<u8>>,
}

pub fn run_annotate(c: &AnnotateConfig, runner: &dyn CommandRunner) -> Result<AnnotateSummary> {
    if c.refdb.is_none() && c.libraries.is_empty() {
        return Err(cli_error("rep annotate needs --refdb and/or --library"));
    }
    fs::create_dir_all(&c.output_dir)?;
    let input = read_input(&c.input)?;
    let db = c.refdb.clone().map(RefDb::new);
    let mut summary = AnnotateSummary {
        sequences: input.seqs.len(),
        total_bp: input.seqs.iter().map(|s| s.len() as u64).sum(),
        ..Default::default()
    };
    let mut commands = Vec::new();
    let mut occupied: HashMap<usize, Occupancy> = HashMap::new();

    // 1. nucleotide hits
    let windows_path = c.output_dir.join("windows.fa");
    let all = (0..input.seqs.len())
        .map(|i| (i, 0, input.seqs[i].len() as u64))
        .collect::<Vec<_>>();
    write_windows(&windows_path, &input, &all, c.window, c.window_overlap)?;
    let mut references: Vec<(String, PathBuf, PathBuf, String)> = Vec::new();
    if let Some(db) = &db {
        references.push((
            "refdb".into(),
            db.nucleotide_library(),
            db.build_dir().join("index"),
            "nucleotides.blastn".into(),
        ));
    }
    for (i, lib) in c.libraries.iter().enumerate() {
        let name = lib
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| format!("library{i}"));
        references.push((
            name,
            lib.clone(),
            c.output_dir.join("index"),
            format!("library{i}.blastn"),
        ));
    }
    let mosaics = match &db {
        Some(db) if db.mosaics().is_file() => crate::refdb::read_mosaics(&db.mosaics())?,
        _ => HashMap::new(),
    };
    let mut nucl = Vec::new();
    for (i, (source, library, index_dir, index_name)) in references.iter().enumerate() {
        fs::create_dir_all(index_dir)?;
        let index = ensure_index(
            index_dir,
            index_name,
            "db",
            &library_stamp(library)?,
            |out| {
                let mut cmd =
                    makeblastdb_nucleotide_command(executable("makeblastdb"), library, out);
                commands.push(command_line(&cmd));
                run_checked(runner, &mut cmd, "makeblastdb").map(drop)
            },
        )?;
        let out = c.output_dir.join(format!("nucleotide_hits.{i}.tsv"));
        let mut cmd = blastn_command(&index, &windows_path, &out, c.threads, c.nucl_evalue);
        commands.push(command_line(&cmd));
        run_checked(runner, &mut cmd, "blastn")?;
        for hit in read_hits(&out)? {
            let mut f = fragment(&hit, source, Evidence::Nucleotide)?;
            // a hit inside an element embedded in a reference of another class
            if source == "refdb" {
                if let Some(m) = crate::refdb::foreign_class(
                    &mosaics,
                    &hit.subject_id,
                    hit.subject_start,
                    hit.subject_end,
                ) {
                    f.class = m.class.clone();
                    f.note = format!("embedded:{}-{}", m.start, m.end);
                }
            }
            nucl.push(f);
        }
    }
    let mut fragments = claim(by_score(nucl), &mut occupied, c.min_fragment);

    // 2. protein hits on what nucleotide hits left
    let protein_library = match (&c.protein_library, &db) {
        (Some(p), _) => Some((
            p.clone(),
            c.output_dir.join("index"),
            "library_proteins.blastx",
        )),
        (None, Some(db)) => Some((
            db.protein_library(),
            db.build_dir().join("index"),
            "proteins.blastx",
        )),
        (None, None) => None,
    };
    if let (true, Some((library, index_dir, index_name))) = (c.proteins, protein_library) {
        let gaps = gaps_of(&input, &occupied, c.protein_min_gap);
        if !gaps.is_empty() {
            let queries = c.output_dir.join("protein_queries.fa");
            write_windows(&queries, &input, &gaps, c.window, c.window_overlap)?;
            fs::create_dir_all(&index_dir)?;
            let index = ensure_index(
                &index_dir,
                index_name,
                "db",
                &library_stamp(&library)?,
                |out| {
                    let mut cmd =
                        makeblastdb_protein_command(executable("makeblastdb"), &library, out);
                    commands.push(command_line(&cmd));
                    run_checked(runner, &mut cmd, "makeblastdb").map(drop)
                },
            )?;
            let out = c.output_dir.join("protein_hits.tsv");
            let mut cmd = blastx_command(&index, &queries, &out, c.threads, c.prot_evalue);
            commands.push(command_line(&cmd));
            run_checked(runner, &mut cmd, "blastx")?;
            let prot = read_hits(&out)?
                .iter()
                .filter(|h| h.bitscore >= c.min_protein_bits)
                .map(|h| fragment(h, "refdb", Evidence::Protein))
                .collect::<Result<Vec<_>>>()?;
            fragments.extend(claim(by_score(prot), &mut occupied, c.min_fragment));
        }
    }

    // 3. tandem and low-complexity runs
    for (scanner, evidence) in [
        (tandem_runs as fn(&[u8]) -> Vec<scan::Run>, Evidence::Tandem),
        (low_complexity_runs, Evidence::LowComplexity),
    ] {
        let mut runs = Vec::new();
        for (seq, start, end) in gaps_of(&input, &occupied, 60) {
            let slice = &input.seqs[seq][start as usize..end as usize];
            for r in scanner(slice) {
                let (class, reference) = match (evidence, r.period) {
                    (Evidence::LowComplexity, _) => ("Low_complexity", "low_complexity".into()),
                    (_, p) if p <= 10 => ("Simple_repeat", format!("({})n", unit(slice, &r))),
                    (_, p) => ("Satellite", format!("tandem_{p}bp")),
                };
                runs.push(Fragment {
                    seq,
                    start: start + r.start as u64,
                    end: start + r.end as u64,
                    strand: '+',
                    reference,
                    source: "scan".into(),
                    class: class.into(),
                    identity: 100.0 * r.identity,
                    evidence,
                    note: if r.period > 0 {
                        format!("period={}", r.period)
                    } else {
                        String::new()
                    },
                    ..Default::default()
                });
            }
        }
        fragments.extend(claim(runs, &mut occupied, c.min_fragment.min(20)));
    }
    fragments.sort_by(|a, b| a.seq.cmp(&b.seq).then(a.start.cmp(&b.start)));
    summary.elements = chain(&mut fragments, c.chain_gap);

    // 4. terminal-repeat pairs in long unannotated stretches
    let mut candidates = Vec::new();
    if c.novel {
        let gaps = gaps_of(&input, &occupied, c.novel_min_gap);
        if !gaps.is_empty() {
            candidates = find_candidates(c, &input, &gaps, runner, &mut commands)?;
            for (k, f) in candidates.iter_mut().enumerate() {
                f.element = summary.elements + k + 1;
            }
        }
    }

    summary.fragments = fragments.len();
    summary.candidate_novel = candidates.len();
    summary.candidate_novel_bp = candidates.iter().map(|f| f.end - f.start).sum();
    for (i, name) in input.names.iter().enumerate() {
        summary.per_sequence.push(SequenceSummary {
            name: name.clone(),
            length: input.seqs[i].len() as u64,
            ..Default::default()
        });
    }
    for f in &fragments {
        let bp = f.end - f.start;
        let class = summary_class(&f.class).to_owned();
        *summary.by_class.entry(class.clone()).or_default() += bp;
        *summary
            .by_evidence
            .entry(f.evidence.as_str().into())
            .or_default() += bp;
        let s = &mut summary.per_sequence[f.seq];
        s.masked_bp += bp;
        *s.by_class.entry(class).or_default() += bp;
        summary.masked_bp += bp;
    }
    summary.hint = single_element_hint(&input, &fragments);
    if let Some(h) = &summary.hint {
        eprintln!("{h}");
    }
    summary.commands = commands;

    let mut all = fragments;
    all.extend(candidates);
    all.sort_by(|a, b| a.seq.cmp(&b.seq).then(a.start.cmp(&b.start)));
    write_table(&c.output_dir.join(TABLE), &input, &all)?;
    write_gff(&c.output_dir.join(GFF), &input, &all)?;
    write_bed(&c.output_dir.join(BED), &input, &all)?;
    write_masked(&c.output_dir.join(MASKED), &input, &all, c.hard_mask)?;
    write_summary_tsv(&c.output_dir.join(SUMMARY_TSV), &summary)?;
    fs::write(
        c.output_dir.join(SUMMARY_JSON),
        serde_json::to_vec_pretty(&summary).map_err(|e| cli_error(e.to_string()))?,
    )?;
    Ok(summary)
}

fn read_input(path: &Path) -> Result<Input> {
    let mut input = Input {
        names: Vec::new(),
        seqs: Vec::new(),
    };
    for record in fasta::Reader::from_file(path)
        .map_err(|e| cli_error(format!("{}: {e}", path.display())))?
        .records()
    {
        let record = record.map_err(|e| cli_error(format!("{}: {e}", path.display())))?;
        input.names.push(record.id().to_owned());
        input.seqs.push(record.seq().to_ascii_uppercase());
    }
    if input.seqs.is_empty() {
        return Err(cli_error(format!("{}: no sequences", path.display())));
    }
    Ok(input)
}

/// Write `[start, end)` stretches of input sequences, cut into overlapping
/// windows, as `q<seq>_<offset>`.
fn write_windows(
    path: &Path,
    input: &Input,
    stretches: &[(usize, u64, u64)],
    window: u64,
    overlap: u64,
) -> Result<()> {
    let mut out = BufWriter::new(fs::File::create(path)?);
    for &(seq, start, end) in stretches {
        for (s, e) in windows(start, end, window, overlap) {
            writeln!(out, ">q{seq}_{s}")?;
            out.write_all(&input.seqs[seq][s as usize..e as usize])?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(())
}

fn windows(start: u64, end: u64, window: u64, overlap: u64) -> Vec<(u64, u64)> {
    let step = window.saturating_sub(overlap).max(1);
    let mut out = Vec::new();
    let mut s = start;
    loop {
        let e = (s + window).min(end);
        out.push((s, e));
        if e >= end {
            break;
        }
        s += step;
    }
    out
}

/// Sequence index and offset of a window ID `q<seq>_<offset>`.
fn window_of(id: &str) -> Result<(usize, u64)> {
    id.strip_prefix('q')
        .and_then(|r| r.split_once('_'))
        .and_then(|(s, o)| Some((s.parse().ok()?, o.parse().ok()?)))
        .ok_or_else(|| cli_error(format!("unexpected query ID in search output: {id}")))
}

fn read_hits(path: &Path) -> Result<Vec<ProteinHit>> {
    if !path.is_file() {
        // the search tools may skip the output file when nothing is found
        return Ok(Vec::new());
    }
    parse_protein_hits(path)
}

fn fragment(hit: &ProteinHit, source: &str, evidence: Evidence) -> Result<Fragment> {
    let (seq, offset) = window_of(&hit.query_id)?;
    // BLASTN reports the minus strand by subject order, BLASTX by query order
    let strand = match evidence {
        Evidence::Nucleotide if hit.subject_start > hit.subject_end => '-',
        Evidence::Nucleotide => '+',
        _ => hit.strand,
    };
    Ok(Fragment {
        seq,
        start: offset + hit.query_start - 1,
        end: offset + hit.query_end,
        strand,
        reference: hit.subject_id.clone(),
        source: source.into(),
        ref_start: hit.subject_start.min(hit.subject_end),
        ref_end: hit.subject_start.max(hit.subject_end),
        ref_length: hit.subject_length,
        class: hit.class.clone().unwrap_or_else(|| "Unknown".into()),
        identity: hit.percent_identity,
        bitscore: hit.bitscore,
        evalue: hit.evalue,
        evidence,
        ..Default::default()
    })
}

/// Hits to classified references first, then best bit score; ties in a
/// fixed order so that runs are reproducible. Unknown references include
/// whole-element entries (Repbase's DM412 is the 412 element, LTRs and all)
/// that would otherwise outscore the classified parts of the same element.
fn by_score(mut hits: Vec<Fragment>) -> Vec<Fragment> {
    let unknown = |f: &Fragment| f.class == "Unknown";
    hits.sort_by(|a, b| {
        unknown(a)
            .cmp(&unknown(b))
            .then(b.bitscore.total_cmp(&a.bitscore))
            .then(a.seq.cmp(&b.seq))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
            .then(a.reference.cmp(&b.reference))
    });
    hits
}

fn gaps_of(
    input: &Input,
    occupied: &HashMap<usize, Occupancy>,
    min: u64,
) -> Vec<(usize, u64, u64)> {
    let empty = Occupancy::default();
    (0..input.seqs.len())
        .flat_map(|i| {
            occupied
                .get(&i)
                .unwrap_or(&empty)
                .gaps(input.seqs[i].len() as u64, min)
                .into_iter()
                .map(move |(s, e)| (i, s, e))
        })
        .collect()
}

/// The repeat unit of a microsatellite run, as written in RepeatMasker
/// (`(CA)n`).
fn unit(slice: &[u8], r: &scan::Run) -> String {
    String::from_utf8_lossy(&slice[r.start..r.start + r.period]).into_owned()
}

/// Self-alignments of the unannotated stretches, read as terminal-repeat pairs.
fn find_candidates(
    c: &AnnotateConfig,
    input: &Input,
    gaps: &[(usize, u64, u64)],
    runner: &dyn CommandRunner,
    commands: &mut Vec<Vec<String>>,
) -> Result<Vec<Fragment>> {
    // long stretches in windows that still hold a 25 kb element whole
    let (window, overlap) = (100_000, 25_000);
    let stretches = c.output_dir.join("unannotated.fa");
    write_windows(&stretches, input, gaps, window, overlap)?;
    let db = c.output_dir.join("unannotated_db").join("db");
    fs::create_dir_all(db.parent().expect("has a parent"))?;
    let mut cmd = makeblastdb_nucleotide_command(executable("makeblastdb"), &stretches, &db);
    commands.push(command_line(&cmd));
    run_checked(runner, &mut cmd, "makeblastdb")?;
    let out = c.output_dir.join("self_hits.tsv");
    let mut cmd = blastn_command(&db, &stretches, &out, c.threads, c.nucl_evalue);
    commands.push(command_line(&cmd));
    run_checked(runner, &mut cmd, "blastn")?;
    let mut by_window: BTreeMap<(usize, u64), Vec<SelfHit>> = BTreeMap::new();
    for h in read_hits(&out)? {
        if h.query_id != h.subject_id {
            continue;
        }
        // the stretch aligned to itself
        if h.query_start == h.subject_start && h.query_end == h.subject_end {
            continue;
        }
        by_window
            .entry(window_of(&h.query_id)?)
            .or_default()
            .push(SelfHit {
                q_start: h.query_start,
                q_end: h.query_end,
                s_start: h.subject_start,
                s_end: h.subject_end,
                identity: h.percent_identity,
            });
    }
    let mut taken: HashMap<usize, Occupancy> = HashMap::new();
    let mut found = Vec::new();
    for ((seq, offset), hits) in by_window {
        for cand in terminal_pairs(&hits) {
            let (start, end) = (offset + cand.start, offset + cand.end);
            let s = &input.seqs[seq];
            let Some(e) = scan::ends(s, cand.kind, start as usize, end as usize) else {
                continue; // no insertion signature
            };
            let (start, end) = (e.start as u64, e.end as u64);
            let occ = taken.entry(seq).or_default();
            if occ.free(start, end) != vec![(start, end)] {
                continue; // found again in an overlapping window
            }
            occ.insert(start, end);
            found.push(Fragment {
                seq,
                start,
                end,
                strand: '.',
                reference: cand.kind.into(),
                source: "scan".into(),
                class: "Unknown".into(),
                identity: cand.identity,
                evidence: Evidence::Structure,
                note: format!(
                    "terminal_repeat={}bp;tsd={}bp{}",
                    cand.repeat_length,
                    e.tsd,
                    if e.tg_ca { ";TG..CA" } else { "" }
                ),
                ..Default::default()
            });
        }
    }
    Ok(found)
}

fn blastn_command(db: &Path, query: &Path, out: &Path, threads: usize, evalue: f64) -> Command {
    let mut c = Command::new(executable("blastn"));
    c.args(["-task", "blastn", "-db"])
        .arg(db)
        .arg("-query")
        .arg(query)
        .args([
            "-outfmt",
            &format!("6 {OUTFMT_FIELDS}"),
            "-evalue",
            &evalue.to_string(),
            // a window can hold hundreds of fragments of related references
            "-max_target_seqs",
            "5000",
            "-num_threads",
            &threads.to_string(),
        ])
        .arg("-out")
        .arg(out);
    c
}

fn blastx_command(db: &Path, query: &Path, out: &Path, threads: usize, evalue: f64) -> Command {
    let mut c = Command::new(executable("blastx"));
    c.arg("-db")
        .arg(db)
        .arg("-query")
        .arg(query)
        .args([
            "-outfmt",
            &format!("6 {OUTFMT_FIELDS}"),
            "-evalue",
            &evalue.to_string(),
            "-max_target_seqs",
            "500",
            "-seg",
            "yes",
            "-num_threads",
            &threads.to_string(),
        ])
        .arg("-out")
        .arg(out);
    c
}

fn command_line(cmd: &Command) -> Vec<String> {
    std::iter::once(cmd.get_program())
        .chain(cmd.get_args())
        .map(|s| s.to_string_lossy().into_owned())
        .collect()
}

/// Top-level class for summaries, with legacy names mapped.
fn summary_class(class: &str) -> &str {
    canonical_class(class)
}

/// A single short sequence mostly covered by one element is better served
/// by `rep classify`.
fn single_element_hint(input: &Input, fragments: &[Fragment]) -> Option<String> {
    if input.seqs.len() != 1 || input.seqs[0].len() > 20_000 {
        return None;
    }
    let len = input.seqs[0].len() as u64;
    let mut by_element: BTreeMap<usize, (u64, &str)> = BTreeMap::new();
    for f in fragments {
        let e = by_element.entry(f.element).or_insert((0, &f.reference));
        e.0 += f.end - f.start;
    }
    let (bp, reference) = by_element.values().max_by_key(|(bp, _)| *bp)?;
    (5 * bp >= 4 * len).then(|| {
        format!(
            "{} looks like a single element ({reference} covers {:.0}%); \
             `rep classify` makes a class call for a consensus",
            input.names[0],
            100.0 * *bp as f64 / len as f64
        )
    })
}

fn write_table(path: &Path, input: &Input, fragments: &[Fragment]) -> Result<()> {
    let mut out = BufWriter::new(fs::File::create(path)?);
    writeln!(
        out,
        "sequence\tstart\tend\tstrand\tclass\treference\tref_start\tref_end\tref_length\tidentity\tbitscore\tevalue\tevidence\tsource\telement\tnote"
    )?;
    for f in fragments {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{:.1}\t{:.2e}\t{}\t{}\t{}\t{}",
            input.names[f.seq],
            f.start + 1,
            f.end,
            f.strand,
            f.class,
            f.reference,
            f.ref_start,
            f.ref_end,
            f.ref_length,
            f.identity,
            f.bitscore,
            f.evalue,
            f.evidence.as_str(),
            f.source,
            f.element,
            f.note
        )?;
    }
    out.flush()?;
    Ok(())
}

/// Sequence Ontology type of a fragment.
fn so_type(f: &Fragment) -> &'static str {
    match f.evidence {
        Evidence::Structure => "repeat_region",
        Evidence::LowComplexity => "low_complexity_region",
        _ => match summary_class(&f.class) {
            "Satellite" => "satellite_DNA",
            "Simple_repeat" => "microsatellite",
            "Low_complexity" => "low_complexity_region",
            c if c == "Unknown" || is_non_te(c) => "repeat_region",
            _ => "transposable_element",
        },
    }
}

/// Escape GFF3 attribute values.
fn gff_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            ';' | '=' | '&' | ',' | '%' | '\t' | '\n' => {
                write!(out, "%{:02X}", ch as u32).expect("writing to a String cannot fail")
            }
            _ => out.push(ch),
        }
    }
    out
}

fn write_gff(path: &Path, input: &Input, fragments: &[Fragment]) -> Result<()> {
    let mut out = BufWriter::new(fs::File::create(path)?);
    writeln!(out, "##gff-version 3")?;
    for (name, seq) in input.names.iter().zip(&input.seqs) {
        writeln!(out, "##sequence-region {name} 1 {}", seq.len())?;
    }
    for (i, f) in fragments.iter().enumerate() {
        let mut attrs = format!(
            "ID=frag{};Name={};class={};element=E{};evidence={}",
            i + 1,
            gff_escape(&f.reference),
            gff_escape(&f.class),
            f.element,
            f.evidence.as_str()
        );
        if f.ref_length > 0 {
            write!(
                attrs,
                ";Target={} {} {};identity={:.1}",
                gff_escape(&f.reference).replace(' ', "%20"),
                f.ref_start,
                f.ref_end,
                f.identity
            )
            .expect("writing to a String cannot fail");
        }
        if !f.note.is_empty() {
            write!(attrs, ";Note={}", gff_escape(&f.note))
                .expect("writing to a String cannot fail");
        }
        let score = if f.bitscore > 0.0 {
            format!("{:.1}", f.bitscore)
        } else {
            ".".into()
        };
        writeln!(
            out,
            "{}\trep\t{}\t{}\t{}\t{score}\t{}\t.\t{attrs}",
            input.names[f.seq],
            so_type(f),
            f.start + 1,
            f.end,
            f.strand
        )?;
    }
    out.flush()?;
    Ok(())
}

fn write_bed(path: &Path, input: &Input, fragments: &[Fragment]) -> Result<()> {
    let mut out = BufWriter::new(fs::File::create(path)?);
    for f in fragments {
        writeln!(
            out,
            "{}\t{}\t{}\t{}|{}\t{}\t{}",
            input.names[f.seq],
            f.start,
            f.end,
            f.class,
            f.reference,
            (f.bitscore.round() as u64).min(1000),
            f.strand
        )?;
    }
    out.flush()?;
    Ok(())
}

/// Soft-masked (lower case) or hard-masked (N) input; candidates are not masked.
fn write_masked(path: &Path, input: &Input, fragments: &[Fragment], hard: bool) -> Result<()> {
    let mut seqs = input.seqs.clone();
    for f in fragments
        .iter()
        .filter(|f| f.evidence != Evidence::Structure)
    {
        for b in &mut seqs[f.seq][f.start as usize..f.end as usize] {
            *b = if hard { b'N' } else { b.to_ascii_lowercase() };
        }
    }
    let mut out = BufWriter::new(fs::File::create(path)?);
    for (name, seq) in input.names.iter().zip(&seqs) {
        writeln!(out, ">{name}")?;
        for line in seq.chunks(60) {
            out.write_all(line)?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(())
}

fn write_summary_tsv(path: &Path, s: &AnnotateSummary) -> Result<()> {
    let mut out = BufWriter::new(fs::File::create(path)?);
    writeln!(out, "sequence\tclass\tbp\tpercent")?;
    let rows = s
        .per_sequence
        .iter()
        .map(|q| (q.name.as_str(), q.length, q.masked_bp, &q.by_class))
        .chain(std::iter::once((
            "total",
            s.total_bp,
            s.masked_bp,
            &s.by_class,
        )));
    for (name, length, masked, classes) in rows {
        let pct = |bp: u64| 100.0 * bp as f64 / length.max(1) as f64;
        for (class, &bp) in classes {
            writeln!(out, "{name}\t{class}\t{bp}\t{:.2}", pct(bp))?;
        }
        writeln!(out, "{name}\tmasked\t{masked}\t{:.2}", pct(masked))?;
    }
    out.flush()?;
    Ok(())
}

fn cli_error(message: impl Into<String>) -> Error {
    Error::new(ErrorKind::GenericCli(message.into()))
}

#[cfg(test)]
mod tests;
