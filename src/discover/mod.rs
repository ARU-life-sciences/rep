//! Cross-species discovery (roadmap M4).
//!
//! Families from many `rep classify --all` runs are pooled and clustered
//! twice with MMseqs2:
//! 1. **Family clusters**: nucleotide clustering at the 80-80 rule
//!    (identity and coverage), collapsing the same element found in many
//!    libraries.
//! 2. **Groups**: family clusters whose ORF proteins cluster together (about
//!    30% identity), linking distant relatives that DNA comparison misses.
//!
//! Labelled members (RepeatClassifier labels, or our own confident calls)
//! are used to propagate a class to Unknown members, with accuracy measured
//! by leave-one-out on the labelled members themselves. Groups with no
//! labelled member are ranked as candidate novel TE classes.

use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fmt::Write as _,
    fs,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    process::Command,
};

use bio::io::fasta;
use serde::Serialize;

mod relatives;
pub use relatives::{relatives_file, write_relatives, RelativesConfig, RelativesSummary};

use crate::{
    classify::{is_non_te, orfs, split_label, top_level},
    command_runner::run_checked,
    CommandRunner, Error, ErrorKind, Result,
};

#[derive(Clone, Debug, Serialize)]
pub struct DiscoverConfig {
    /// Root directory searched for `rep classify --all` outputs
    /// (`classification.tsv` with `queries.fa` beside it).
    pub results: PathBuf,
    /// Optional TSV with `scientific_name` and `class` columns.
    pub taxa: Option<PathBuf>,
    pub output_dir: PathBuf,
    pub threads: usize,
    pub family_identity: f64,
    pub family_coverage: f64,
    pub protein_identity: f64,
    pub protein_coverage: f64,
    pub min_orf_aa: usize,
    /// ORFs per family used for protein clustering (longest first).
    pub orfs_per_family: usize,
    /// Minimum fraction of labelled members agreeing to propagate a label.
    pub min_purity: f64,
}

impl DiscoverConfig {
    pub fn new(results: PathBuf, output_dir: PathBuf) -> Self {
        Self {
            results,
            taxa: None,
            output_dir,
            threads: 1,
            family_identity: 0.8,
            family_coverage: 0.8,
            protein_identity: 0.3,
            protein_coverage: 0.5,
            min_orf_aa: 100,
            orfs_per_family: 3,
            min_purity: 0.75,
        }
    }
}

/// One pooled family.
#[derive(Clone, Debug, Default)]
struct Member {
    sample: String,
    species: String,
    taxon: String,
    family_id: String,
    length: usize,
    original_label: String,
    disposition: String,
    label: String,
    /// Label trusted as a reference for propagation, if any.
    reference: Option<String>,
    hallmarks: Vec<&'static str>,
    genome_bp: u64,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct DiscoverSummary {
    pub samples: usize,
    pub members: usize,
    pub labelled_members: usize,
    pub family_clusters: usize,
    pub groups: usize,
    pub unlabelled_members: usize,
    pub propagated_from_family_cluster: usize,
    pub propagated_from_group: usize,
    /// Members whose propagated label differs in class from their own
    /// (`relabels.tsv`, for `rep relabel --overrides`).
    pub relabels: usize,
    /// Leave-one-out accuracy of propagation on labelled members, by level.
    pub loo_family_cluster: Accuracy,
    pub loo_group: Accuracy,
    pub novel_groups: usize,
    /// Novel groups with members from at least two species.
    pub novel_groups_multi_species: usize,
    /// Protein clusters not used for linking because their labelled
    /// members span more than one top-level class.
    pub mixed_protein_clusters: usize,
    /// Protein links refused because they would join components labelled
    /// with different classes.
    pub blocked_links: usize,
    pub commands: Vec<Vec<String>>,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct Accuracy {
    pub predicted: usize,
    pub correct: usize,
    pub accuracy: f64,
}

pub fn run_discover(c: &DiscoverConfig, runner: &dyn CommandRunner) -> Result<DiscoverSummary> {
    fs::create_dir_all(&c.output_dir)?;
    let taxa = match &c.taxa {
        Some(p) => read_taxa(p, "class")?,
        None => HashMap::new(),
    };
    let (members, samples) = pool(c, &taxa)?;
    eprintln!("Pooled {} families from {samples} samples", members.len());
    let mut summary = DiscoverSummary {
        samples,
        members: members.len(),
        labelled_members: members.iter().filter(|m| m.reference.is_some()).count(),
        ..Default::default()
    };

    let mut commands = Vec::new();
    let nucl_clusters = mmseqs_cluster(
        c,
        runner,
        "linclust",
        &c.output_dir.join("pooled.fa"),
        &c.output_dir.join("families"),
        c.family_identity,
        c.family_coverage,
        &mut commands,
    )?;
    let prot_path = c.output_dir.join("pooled_orfs.faa");
    let prot_clusters = if fs::metadata(&prot_path)?.len() > 0 {
        mmseqs_cluster(
            c,
            runner,
            "cluster",
            &prot_path,
            &c.output_dir.join("proteins"),
            c.protein_identity,
            c.protein_coverage,
            &mut commands,
        )?
    } else {
        Vec::new()
    };
    summary.commands = commands;

    // family clusters: member index -> cluster representative index
    let index = |id: &str| {
        id.trim_start_matches('m')
            .split('_')
            .next()
            .and_then(|n| n.parse::<usize>().ok())
    };
    let mut family_of = (0..members.len()).collect::<Vec<_>>();
    for (rep, member) in &nucl_clusters {
        if let (Some(r), Some(m)) = (index(rep), index(member)) {
            family_of[m] = r;
        }
    }
    // groups: union family clusters that share a protein cluster
    // each component records the top-level classes of its labelled members
    let mut class_bits: BTreeMap<&str, u32> = BTreeMap::new();
    let masks = members
        .iter()
        .map(|m| match m.reference.as_deref().map(top_level) {
            Some(class) => {
                let next = class_bits.len() as u32;
                1u32 << (*class_bits.entry(class).or_insert(next.min(31)))
            }
            None => 0,
        })
        .collect::<Vec<_>>();
    let mut uf = UnionFind::new(masks);
    for (m, &f) in family_of.iter().enumerate() {
        uf.union(m, f);
    }
    // a protein cluster whose labelled members disagree on the class is a
    // bridge (captured gene fragments, chimeric families), not evidence of
    // relatedness, so it is not used for linking
    let mut protein_members: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (rep, member) in &prot_clusters {
        if let Some(m) = index(member) {
            protein_members.entry(rep.as_str()).or_default().push(m);
        }
    }
    for ids in protein_members.values() {
        let classes = ids
            .iter()
            .filter_map(|&i| members[i].reference.as_deref().map(top_level))
            .collect::<BTreeSet<_>>();
        if classes.len() > 1 {
            summary.mixed_protein_clusters += 1;
            continue;
        }
        for &m in &ids[1..] {
            // cannot-link: never join components labelled with different classes,
            // which chains of unlabelled protein clusters would otherwise do
            if !uf.union_if_compatible(ids[0], m) {
                summary.blocked_links += 1;
            }
        }
    }
    let group_of = (0..members.len()).map(|m| uf.find(m)).collect::<Vec<_>>();
    summary.family_clusters = family_of.iter().collect::<BTreeSet<_>>().len();
    summary.groups = group_of.iter().collect::<BTreeSet<_>>().len();

    // label propagation with leave-one-out accuracy
    let by_family = invert(&family_of);
    let by_group = invert(&group_of);
    let mut propagated: Vec<Option<(String, &'static str, f64, usize)>> = vec![None; members.len()];
    for (m, member) in members.iter().enumerate() {
        let family = predict(&members, &by_family[&family_of[m]], m, c.min_purity);
        let group = predict(&members, &by_group[&group_of[m]], m, c.min_purity);
        match &member.reference {
            Some(truth) => {
                for (prediction, acc) in [
                    (&family, &mut summary.loo_family_cluster),
                    (&group, &mut summary.loo_group),
                ] {
                    if let Some((label, _, _)) = prediction {
                        acc.predicted += 1;
                        acc.correct += usize::from(top_level(label) == top_level(truth));
                    }
                }
            }
            None => {
                if let Some((label, purity, n)) = family {
                    propagated[m] = Some((label, "family_cluster", purity, n));
                    summary.propagated_from_family_cluster += 1;
                } else if let Some((label, purity, n)) = group {
                    propagated[m] = Some((label, "group", purity, n));
                    summary.propagated_from_group += 1;
                }
            }
        }
    }
    for acc in [&mut summary.loo_family_cluster, &mut summary.loo_group] {
        acc.accuracy = acc.correct as f64 / acc.predicted.max(1) as f64;
    }
    summary.unlabelled_members = members.iter().filter(|m| m.reference.is_none()).count();

    write_members(c, &members, &family_of, &group_of, &propagated)?;
    summary.relabels = write_relabels(c, &members, &propagated)?;
    (summary.novel_groups, summary.novel_groups_multi_species) =
        write_groups(c, &members, &by_group, &family_of)?;
    fs::write(
        c.output_dir.join("summary.json"),
        serde_json::to_vec_pretty(&summary).map_err(|e| cli_error(e.to_string()))?,
    )?;
    Ok(summary)
}

/// Majority top-level label among the other labelled members, if pure enough.
fn predict(
    members: &[Member],
    cluster: &[usize],
    exclude: usize,
    min_purity: f64,
) -> Option<(String, f64, usize)> {
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    let mut full: BTreeMap<&str, BTreeMap<&str, usize>> = BTreeMap::new();
    let mut n = 0;
    for &i in cluster.iter().filter(|&&i| i != exclude) {
        if let Some(label) = &members[i].reference {
            *counts.entry(top_level(label)).or_default() += 1;
            *full
                .entry(top_level(label))
                .or_default()
                .entry(label)
                .or_default() += 1;
            n += 1;
        }
    }
    let (best, count) = counts
        .into_iter()
        .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(a.0)))?;
    let purity = count as f64 / n as f64;
    if purity < min_purity {
        return None;
    }
    // the most common full label within the winning top-level class
    let label = full[best]
        .iter()
        .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))
        .map(|(l, _)| (*l).to_owned())?;
    Some((label, purity, n))
}

fn invert(assignment: &[usize]) -> HashMap<usize, Vec<usize>> {
    let mut map: HashMap<usize, Vec<usize>> = HashMap::new();
    for (m, &c) in assignment.iter().enumerate() {
        map.entry(c).or_default().push(m);
    }
    map
}

/// Pool every family of every classify output into one nucleotide FASTA and
/// one ORF protein FASTA with compact IDs (`m<index>`, `m<index>_o<k>`).
fn pool(c: &DiscoverConfig, taxa: &HashMap<String, String>) -> Result<(Vec<Member>, usize)> {
    let dirs = classify_dirs(&c.results);
    let mut nucl = BufWriter::new(fs::File::create(c.output_dir.join("pooled.fa"))?);
    let mut prot = BufWriter::new(fs::File::create(c.output_dir.join("pooled_orfs.faa"))?);
    let mut members = Vec::new();
    for dir in &dirs {
        let sample = dir
            .strip_prefix(&c.results)
            .unwrap_or(dir)
            .to_string_lossy()
            .into_owned();
        let species = sample.split('/').next().unwrap_or(&sample).to_owned();
        let taxon = taxa.get(&species).cloned().unwrap_or_default();
        let rows = read_classification(&dir.join("classification.tsv"))?;
        for record in fasta::Reader::from_file(dir.join("queries.fa"))?.records() {
            let record = record?;
            let Some(row) = rows.get(record.id()) else {
                continue;
            };
            let idx = members.len();
            let seq = record.seq().to_ascii_uppercase();
            writeln!(nucl, ">m{idx}")?;
            nucl.write_all(&seq)?;
            nucl.write_all(b"\n")?;
            for (k, orf) in orfs(&seq, c.min_orf_aa)
                .iter()
                .take(c.orfs_per_family)
                .enumerate()
            {
                writeln!(prot, ">m{idx}_o{k}\n{}", orf.protein)?;
            }
            members.push(member_from_row(row, &sample, &species, &taxon, seq.len()));
        }
    }
    nucl.flush()?;
    prot.flush()?;
    Ok((members, dirs.len()))
}

/// `rep classify` output directories under `results` (`classification.tsv` with `queries.fa`).
fn classify_dirs(results: &Path) -> Vec<PathBuf> {
    let mut dirs = walkdir::WalkDir::new(results)
        .follow_links(true)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name() == "classification.tsv")
        .filter_map(|e| e.path().parent().map(Path::to_path_buf))
        .filter(|d| d.join("queries.fa").is_file())
        .collect::<Vec<_>>();
    dirs.sort();
    dirs
}

fn read_classification(path: &Path) -> Result<HashMap<String, HashMap<String, String>>> {
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    let headers = reader.headers()?.clone();
    let mut rows = HashMap::new();
    for row in reader.records() {
        let row = row?;
        let map = headers
            .iter()
            .zip(row.iter())
            .map(|(h, v)| (h.to_owned(), v.to_owned()))
            .collect::<HashMap<_, _>>();
        if let Some(id) = map.get("family_id").cloned() {
            rows.insert(id, map);
        }
    }
    Ok(rows)
}

fn member_from_row(
    row: &HashMap<String, String>,
    sample: &str,
    species: &str,
    taxon: &str,
    length: usize,
) -> Member {
    let get = |k: &str| row.get(k).cloned().unwrap_or_default();
    let original = get("original_label");
    let (disposition, label, confidence) = (get("disposition"), get("label"), get("confidence"));
    // a reference label: our own confident TE call (the learned classifier's
    // decision, or a cascade call above low confidence), else RepeatClassifier's
    // TE label where we did not contradict its class; never a label from relatives
    let usable = |l: &str| !l.is_empty() && !l.starts_with("Unknown") && !is_non_te(l);
    let decided_by = get("decided_by");
    let reference = match decided_by.as_str() {
        "relative" => None,
        "model" if usable(&label) => Some(label.clone()),
        _ if matches!(disposition.as_str(), "known_te" | "chimeric")
            && confidence != "low"
            && usable(&label) =>
        {
            Some(label.clone())
        }
        _ if usable(&original) && top_level(&label) == top_level(&original) => {
            Some(original.clone())
        }
        _ => None,
    };
    let mut hallmarks = Vec::new();
    for (key, name) in [
        ("ltr_length", "LTR"),
        ("tir_length", "TIR"),
        ("tsd_length", "TSD"),
        ("domain_novelty", "novel_domains"),
        ("domain_architecture", "domains"),
    ] {
        if !get(key).is_empty() {
            hallmarks.push(name);
        }
    }
    if get("longest_orf_aa").parse::<usize>().unwrap_or(0) >= 300 {
        hallmarks.push("ORF300");
    }
    Member {
        sample: sample.to_owned(),
        species: species.to_owned(),
        taxon: taxon.to_owned(),
        family_id: get("family_id"),
        length,
        original_label: original,
        disposition,
        label,
        reference,
        hallmarks,
        genome_bp: get("genome_bp").parse().unwrap_or(0),
    }
}

/// `scientific_name` → `class` (falling back to `phylum`), keyed by the
/// lower-case, underscore-separated species name used for directories.
/// Species (lower case, `_` for spaces) to their taxon at `rank` (a column
/// such as `class` or `order`), falling back to the phylum when it is empty.
fn read_taxa(path: &Path, rank: &str) -> Result<HashMap<String, String>> {
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    let headers = reader.headers()?.clone();
    let col = |n: &str| headers.iter().position(|h| h.trim_matches('"') == n);
    let (Some(name), Some(class)) = (col("scientific_name").or_else(|| col("species")), col(rank))
    else {
        return Err(cli_error(format!(
            "{} needs scientific_name (or species) and {rank} columns",
            path.display()
        )));
    };
    let phylum = col("phylum");
    let mut taxa = HashMap::new();
    for row in reader.records() {
        let row = row?;
        let clean = |s: &str| s.trim_matches('"').to_owned();
        let mut value = clean(&row[class]);
        if value.is_empty() {
            value = phylum.map(|p| clean(&row[p])).unwrap_or_default();
        }
        taxa.insert(clean(&row[name]).to_lowercase().replace(' ', "_"), value);
    }
    Ok(taxa)
}

#[allow(clippy::too_many_arguments)]
fn mmseqs_cluster(
    c: &DiscoverConfig,
    runner: &dyn CommandRunner,
    mode: &str,
    input: &Path,
    prefix: &Path,
    identity: f64,
    coverage: f64,
    commands: &mut Vec<Vec<String>>,
) -> Result<Vec<(String, String)>> {
    let tmp = c.output_dir.join(format!("mmseqs_tmp_{mode}"));
    let mut cmd = Command::new(crate::copy_recovery::executable("mmseqs"));
    cmd.arg(format!("easy-{mode}"))
        .arg(input)
        .arg(prefix)
        .arg(&tmp)
        .args([
            "--min-seq-id",
            &identity.to_string(),
            "-c",
            &coverage.to_string(),
        ])
        .args([
            "--cov-mode",
            "0",
            "-v",
            "1",
            "--threads",
            &c.threads.to_string(),
        ]);
    commands.push(
        std::iter::once(cmd.get_program())
            .chain(cmd.get_args())
            .map(|s| s.to_string_lossy().into_owned())
            .collect(),
    );
    run_checked(runner, &mut cmd, &format!("mmseqs easy-{mode}"))?;
    fs::remove_dir_all(&tmp).ok();
    let path = PathBuf::from(format!("{}_cluster.tsv", prefix.display()));
    let text =
        fs::read_to_string(&path).map_err(|e| cli_error(format!("{}: {e}", path.display())))?;
    Ok(text
        .lines()
        .filter_map(|l| l.split_once('\t'))
        .map(|(a, b)| (a.to_owned(), b.to_owned()))
        .collect())
}

fn write_members(
    c: &DiscoverConfig,
    members: &[Member],
    family_of: &[usize],
    group_of: &[usize],
    propagated: &[Option<(String, &'static str, f64, usize)>],
) -> Result<()> {
    let mut out = BufWriter::new(fs::File::create(c.output_dir.join("members.tsv"))?);
    writeln!(
        out,
        "member\tsample\ttaxon\tfamily_id\tlength\toriginal_label\tdisposition\tlabel\treference_label\tfamily_cluster\tgroup\tpropagated_label\tpropagated_from\tpurity\tlabelled_neighbours\thallmarks\tgenome_bp"
    )?;
    for (i, m) in members.iter().enumerate() {
        let (pl, from, purity, n) = match &propagated[i] {
            Some((l, f, p, n)) => (l.as_str(), *f, format!("{p:.2}"), n.to_string()),
            None => ("", "", String::new(), String::new()),
        };
        writeln!(
            out,
            "m{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tm{}\tm{}\t{pl}\t{from}\t{purity}\t{n}\t{}\t{}",
            m.sample,
            m.taxon,
            m.family_id,
            m.length,
            m.original_label,
            m.disposition,
            m.label,
            m.reference.as_deref().unwrap_or(""),
            family_of[i],
            group_of[i],
            m.hallmarks.join(","),
            m.genome_bp
        )?;
    }
    out.flush()?;
    Ok(())
}

/// Members whose propagated label differs in class from their own label,
/// for `rep relabel --overrides`. Propagation reaches only members without a
/// reference label, so these are Unknowns and labels from relatives or at
/// low confidence, corrected by the confident labels of their cluster.
fn write_relabels(
    c: &DiscoverConfig,
    members: &[Member],
    propagated: &[Option<(String, &'static str, f64, usize)>],
) -> Result<usize> {
    let mut out = BufWriter::new(fs::File::create(c.output_dir.join(RELABELS))?);
    writeln!(
        out,
        "sample\tfamily_id\tlabel\tnew_label\tsource\tpurity\tlabelled_neighbours"
    )?;
    let mut n = 0;
    for (m, p) in members.iter().zip(propagated) {
        let Some((label, from, purity, k)) = p else {
            continue;
        };
        if top_level(label) == top_level(&m.label) {
            continue;
        }
        writeln!(
            out,
            "{}\t{}\t{}\t{label}\t{from}\t{purity:.2}\t{k}",
            m.sample, m.family_id, m.label
        )?;
        n += 1;
    }
    out.flush()?;
    Ok(n)
}

/// Relabels written by `rep discover`.
pub const RELABELS: &str = "relabels.tsv";

/// Write group summaries, ranking groups without any labelled member as
/// candidate novel classes. Returns the number of such groups, and of those
/// spanning at least two species.
fn write_groups(
    c: &DiscoverConfig,
    members: &[Member],
    by_group: &HashMap<usize, Vec<usize>>,
    family_of: &[usize],
) -> Result<(usize, usize)> {
    let mut rows = Vec::new();
    for (&g, ids) in by_group {
        let species = ids
            .iter()
            .map(|&i| members[i].species.as_str())
            .collect::<BTreeSet<_>>();
        let taxa = ids
            .iter()
            .map(|&i| members[i].taxon.as_str())
            .filter(|t| !t.is_empty())
            .collect::<BTreeSet<_>>();
        let clusters = ids.iter().map(|&i| family_of[i]).collect::<BTreeSet<_>>();
        let mut labels: BTreeMap<&str, usize> = BTreeMap::new();
        for &i in ids {
            if let Some(l) = &members[i].reference {
                *labels.entry(top_level(l)).or_default() += 1;
            }
        }
        let labelled: usize = labels.values().sum();
        // members already explained as simple, satellite or host sequence
        let non_te = ids
            .iter()
            .filter(|&&i| {
                let m = &members[i];
                is_non_te(&m.label) || is_non_te(&m.original_label)
            })
            .count();
        let mut hallmarks: BTreeMap<&str, usize> = BTreeMap::new();
        for &i in ids {
            for h in &members[i].hallmarks {
                *hallmarks.entry(h).or_default() += 1;
            }
        }
        let genome_bp: u64 = ids.iter().map(|&i| members[i].genome_bp).sum();
        let frac = |h: &str| hallmarks.get(h).copied().unwrap_or(0) as f64 / ids.len() as f64;
        // a novel class has no TE-labelled member and is not mostly non-TE
        let novel = labelled == 0 && non_te * 2 < ids.len();
        // breadth across species dominates, then abundance and TE hallmarks;
        // only meaningful for novel groups
        let score = 2.0 * (1.0 + species.len() as f64).ln()
            + 0.5 * (1.0 + ids.len() as f64).ln()
            + frac("TIR").max(frac("LTR")) * 2.0
            + frac("TSD") * 2.0
            + if hallmarks.contains_key("ORF300") {
                1.0
            } else {
                0.0
            }
            + if hallmarks.contains_key("novel_domains") {
                2.0
            } else {
                0.0
            }
            + (1.0 + genome_bp as f64).log10() / 2.0;
        let mut line = String::new();
        write!(
            line,
            "m{g}\t{}\t{}\t{}\t{}\t{}\t{labelled}\t{}\t{non_te}\t{}\t{genome_bp}\t{novel}\t{score:.3}\t{}",
            ids.len(),
            species.len(),
            clusters.len(),
            taxa.len(),
            taxa.iter().copied().collect::<Vec<_>>().join(","),
            labels
                .iter()
                .map(|(l, n)| format!("{l}:{n}"))
                .collect::<Vec<_>>()
                .join(","),
            hallmarks
                .iter()
                .map(|(h, n)| format!("{h}:{n}"))
                .collect::<Vec<_>>()
                .join(","),
            ids.iter()
                .take(5)
                .map(|&i| format!("{}|{}", members[i].species, members[i].family_id))
                .collect::<Vec<_>>()
                .join(";"),
        )
        .expect("writing to a String cannot fail");
        rows.push((novel, score, line, species.len()));
    }
    rows.sort_by(|a, b| b.0.cmp(&a.0).then(b.1.total_cmp(&a.1)));
    let mut out = BufWriter::new(fs::File::create(c.output_dir.join("groups.tsv"))?);
    writeln!(
        out,
        "group\tmembers\tspecies\tfamily_clusters\ttaxa\ttaxon_list\tlabelled\tlabels\tnon_te\thallmarks\tgenome_bp\tnovel\tnovelty_score\texamples"
    )?;
    for (_, _, line, _) in &rows {
        writeln!(out, "{line}")?;
    }
    out.flush()?;
    Ok((
        rows.iter().filter(|r| r.0).count(),
        rows.iter().filter(|r| r.0 && r.3 >= 2).count(),
    ))
}

/// Union-find whose components carry a bit mask of labelled classes.
struct UnionFind {
    parent: Vec<usize>,
    mask: Vec<u32>,
}

impl UnionFind {
    fn new(mask: Vec<u32>) -> Self {
        Self {
            parent: (0..mask.len()).collect(),
            mask,
        }
    }
    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            let (root, child) = (ra.min(rb), ra.max(rb));
            self.parent[child] = root;
            self.mask[root] |= self.mask[child];
        }
    }
    /// Join unless both components are labelled and their classes differ.
    fn union_if_compatible(&mut self, a: usize, b: usize) -> bool {
        let (ra, rb) = (self.find(a), self.find(b));
        let (ma, mb) = (self.mask[ra], self.mask[rb]);
        if ra == rb || ma == 0 || mb == 0 || ma == mb {
            self.union(ra, rb);
            true
        } else {
            false
        }
    }
}

fn cli_error(message: impl Into<String>) -> Error {
    Error::new(ErrorKind::GenericCli(message.into()))
}

/// Split `name#label` for callers that only have a family ID.
pub fn family_name(id: &str) -> String {
    split_label(id).0
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::process::ExitStatusExt;
    use std::process::{ExitStatus, Output};

    /// Fakes MMseqs2: clusters pooled members by a fixed assignment.
    struct FakeMmseqs {
        nucl: Vec<(usize, usize)>,
        prot: Vec<(String, String)>,
    }

    impl CommandRunner for FakeMmseqs {
        fn run(&self, cmd: &mut Command) -> Result<Output> {
            let args = cmd
                .get_args()
                .map(|a| a.to_string_lossy().into_owned())
                .collect::<Vec<_>>();
            let prefix = &args[2];
            let body = if args[0] == "easy-linclust" {
                self.nucl
                    .iter()
                    .map(|(r, m)| format!("m{r}\tm{m}\n"))
                    .collect::<String>()
            } else {
                self.prot
                    .iter()
                    .map(|(r, m)| format!("{r}\t{m}\n"))
                    .collect::<String>()
            };
            fs::write(format!("{prefix}_cluster.tsv"), body)?;
            Ok(Output {
                status: ExitStatus::from_raw(0),
                stdout: vec![],
                stderr: vec![],
            })
        }
    }

    fn sample(dir: &Path, name: &str, rows: &[(&str, &str, &str, &str, &str)]) {
        let d = dir.join(name);
        fs::create_dir_all(&d).unwrap();
        let mut tsv = String::from("family_id\toriginal_label\tdisposition\tlabel\tconfidence\ttir_length\tlongest_orf_aa\tgenome_bp\n");
        let mut fa = String::new();
        for (id, orig, disp, label, conf) in rows {
            tsv.push_str(&format!(
                "{id}\t{orig}\t{disp}\t{label}\t{conf}\t\t0\t1000\n"
            ));
            // 120 codons without a stop: one ORF per family for protein linking
            fa.push_str(&format!(">{id}\n{}\n", "GCT".repeat(120)));
        }
        fs::write(d.join("classification.tsv"), tsv).unwrap();
        fs::write(d.join("queries.fa"), fa).unwrap();
    }

    #[test]
    fn references_are_confident_rep_calls_or_uncontradicted_repeatmodeler_labels() {
        let reference = |orig: &str, disp: &str, label: &str, conf: &str, by: &str| {
            let row = [
                ("original_label", orig),
                ("disposition", disp),
                ("label", label),
                ("confidence", conf),
                ("decided_by", by),
            ]
            .into_iter()
            .map(|(k, v)| (k.to_owned(), v.to_owned()))
            .collect::<HashMap<_, _>>();
            member_from_row(&row, "s", "s", "t", 100).reference
        };
        let some = |l: &str| Some(l.to_owned());
        assert_eq!(
            reference(
                "Unknown",
                "interspersed_repeat",
                "DNA/hAT",
                "medium",
                "relative"
            ),
            None
        );
        assert_eq!(
            reference("Unknown", "unresolved", "LINE/L2", "medium", "model"),
            some("LINE/L2")
        );
        assert_eq!(
            reference("SINE/tRNA", "known_te", "LINE/RTE", "high", "cascade"),
            some("LINE/RTE")
        );
        // RepeatModeler's label stands where we agree on the class or found nothing
        assert_eq!(
            reference(
                "LTR/Gypsy",
                "interspersed_repeat",
                "LTR/Gypsy",
                "low",
                "cascade"
            ),
            some("LTR/Gypsy")
        );
        // but not where we called another class, even at low confidence
        assert_eq!(
            reference("SINE/tRNA", "known_te", "DNA/hAT", "low", "cascade"),
            None
        );
    }

    #[test]
    fn propagates_labels_and_ranks_novel_groups() {
        let dir = tempfile::tempdir().unwrap();
        let results = dir.path().join("results");
        // species a: two labelled LINEs and an Unknown; species b: Unknowns
        sample(
            &results,
            "sp_a/GCA_1",
            &[
                ("f1#LINE/L2", "LINE/L2", "known_te", "LINE/L2", "high"),
                ("f2#LINE/L2", "LINE/L2", "known_te", "LINE/L2", "high"),
                ("f3#Unknown", "Unknown", "unresolved", "Unknown", "low"),
            ],
        );
        sample(
            &results,
            "sp_b/GCA_2",
            &[
                ("g1#Unknown", "Unknown", "unresolved", "Unknown", "low"),
                (
                    "g2#Unknown",
                    "Unknown",
                    "novel_candidate",
                    "DNA/Unknown",
                    "medium",
                ),
                (
                    "g3#Unknown",
                    "Unknown",
                    "interspersed_repeat",
                    "Unknown",
                    "medium",
                ),
                (
                    "h1#Simple_repeat",
                    "Simple_repeat",
                    "simple_repeat",
                    "Simple_repeat",
                    "medium",
                ),
                ("h2#Unknown", "Unknown", "unresolved", "Unknown", "low"),
            ],
        );
        // members: 0..2 sp_a, 3..7 sp_b. Family clusters: {0,1,2,3}, {4}, {5}
        // and {6,7}, half simple repeat (not novel). Proteins link 4 and 5.
        let runner = FakeMmseqs {
            nucl: vec![
                (0, 0),
                (0, 1),
                (0, 2),
                (0, 3),
                (4, 4),
                (5, 5),
                (6, 6),
                (6, 7),
            ],
            prot: vec![
                ("m4_o0".into(), "m4_o0".into()),
                ("m4_o0".into(), "m5_o0".into()),
            ],
        };
        let config = DiscoverConfig::new(results, dir.path().join("out"));
        let s = run_discover(&config, &runner).unwrap();
        assert_eq!((s.samples, s.members, s.labelled_members), (2, 8, 2));
        assert_eq!((s.family_clusters, s.groups), (4, 3));
        assert_eq!(
            s.propagated_from_family_cluster, 2,
            "f3 and g1 inherit LINE"
        );
        // each labelled LINE is predicted from the other one: 2/2 correct
        assert_eq!(
            (s.loo_family_cluster.predicted, s.loo_family_cluster.correct),
            (2, 2)
        );
        assert_eq!((s.novel_groups, s.novel_groups_multi_species), (1, 0));
        assert_eq!(s.relabels, 2);
        let relabels = fs::read_to_string(dir.path().join("out").join(RELABELS)).unwrap();
        assert!(relabels.contains("sp_b/GCA_2\tg1#Unknown\tUnknown\tLINE/L2\tfamily_cluster\t"));
        let members = fs::read_to_string(dir.path().join("out/members.tsv")).unwrap();
        assert!(members.lines().any(
            |l| l.contains("g1#Unknown") && l.contains("\tLINE/L2\tfamily_cluster\t1.00\t2\t")
        ));
        let groups = fs::read_to_string(dir.path().join("out/groups.tsv")).unwrap();
        let first = groups.lines().nth(1).unwrap();
        assert!(first.starts_with("m4\t2\t1\t2\t"), "{first}");
        assert!(first.contains("\ttrue\t"));
        assert!(groups
            .lines()
            .any(|l| l.starts_with("m6\t2\t") && l.contains("\t1\t\t") && l.contains("\tfalse\t")));
    }

    #[test]
    fn purity_threshold_blocks_mixed_clusters() {
        let members = ["LINE/L2", "LTR/Gypsy", "LINE/L1"]
            .iter()
            .map(|l| Member {
                reference: Some((*l).into()),
                ..Default::default()
            })
            .chain(std::iter::once(Member::default()))
            .collect::<Vec<_>>();
        let all = [0, 1, 2, 3];
        assert_eq!(
            predict(&members, &all, 3, 0.75),
            None,
            "LINE 2/3 is below 0.75"
        );
        let (label, purity, n) = predict(&members, &all, 3, 0.6).unwrap();
        assert_eq!((top_level(&label), n), ("LINE", 3));
        assert!((purity - 2.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn protein_clusters_spanning_classes_do_not_link() {
        let dir = tempfile::tempdir().unwrap();
        let results = dir.path().join("results");
        sample(
            &results,
            "sp_a/GCA_1",
            &[
                ("f1#LINE/L2", "LINE/L2", "known_te", "LINE/L2", "high"),
                ("f2#LTR/Gypsy", "LTR/Gypsy", "known_te", "LTR/Gypsy", "high"),
                ("f3#Unknown", "Unknown", "unresolved", "Unknown", "low"),
            ],
        );
        let runner = FakeMmseqs {
            nucl: vec![(0, 0), (1, 1), (2, 2)],
            // one protein cluster contains a LINE, an LTR and the Unknown
            prot: vec![
                ("m0_o0".into(), "m0_o0".into()),
                ("m0_o0".into(), "m1_o0".into()),
                ("m0_o0".into(), "m2_o0".into()),
            ],
        };
        let s = run_discover(
            &DiscoverConfig::new(results, dir.path().join("out")),
            &runner,
        )
        .unwrap();
        assert_eq!((s.groups, s.mixed_protein_clusters), (3, 1));
        assert_eq!(s.propagated_from_group, 0);
    }

    #[test]
    fn unlabelled_bridges_cannot_join_differently_labelled_components() {
        let dir = tempfile::tempdir().unwrap();
        let results = dir.path().join("results");
        sample(
            &results,
            "sp_a/GCA_1",
            &[
                ("f1#LINE/L2", "LINE/L2", "known_te", "LINE/L2", "high"),
                ("f2#Unknown", "Unknown", "unresolved", "Unknown", "low"),
                ("f3#LTR/Gypsy", "LTR/Gypsy", "known_te", "LTR/Gypsy", "high"),
            ],
        );
        // two unlabelled-looking protein clusters chain LINE - Unknown - LTR
        let runner = FakeMmseqs {
            nucl: vec![(0, 0), (1, 1), (2, 2)],
            prot: vec![
                ("m0_o0".into(), "m0_o0".into()),
                ("m0_o0".into(), "m1_o0".into()),
                ("m2_o0".into(), "m2_o0".into()),
                ("m2_o0".into(), "m1_o0".into()),
            ],
        };
        let s = run_discover(
            &DiscoverConfig::new(results, dir.path().join("out")),
            &runner,
        )
        .unwrap();
        // the Unknown joins the first component it meets; the LTR stays separate
        assert_eq!((s.groups, s.blocked_links), (2, 1));
    }
}
