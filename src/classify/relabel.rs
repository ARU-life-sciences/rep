//! `rep relabel`: carry `rep classify` labels into an existing RepeatMasker
//! annotation, and measure what changed in base pairs.
//!
//! RepeatMasker names the library family on every `.out` line, so the
//! relabelled library does not have to be masked again: the class/family
//! column is rewritten from `classification.tsv`. Base pairs are counted as
//! the union of each class's intervals, so overlapping hits count once.

use std::{
    collections::{BTreeMap, HashMap},
    fs,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use serde::Serialize;

use crate::{Error, ErrorKind, Result};

#[derive(Clone, Debug)]
pub struct RelabelConfig {
    /// A `rep classify` output directory (`classification.tsv`).
    pub classify_dir: PathBuf,
    pub rm_out: PathBuf,
    /// `.sizes` or `.fai`: sequence names and lengths, for percentages.
    pub genome_sizes: Option<PathBuf>,
    /// `rep discover` relabels (`relabels.tsv`), applied over the classify
    /// labels for rows of `sample` (all rows without one).
    pub overrides: Option<PathBuf>,
    pub sample: Option<String>,
    pub output_dir: PathBuf,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct RelabelSummary {
    pub hits: usize,
    pub relabelled_hits: usize,
    pub genome_bp: Option<u64>,
    pub masked_bp: u64,
    /// Class (the label's first component) to base pairs, before and after.
    pub bp_before: BTreeMap<String, u64>,
    pub bp_after: BTreeMap<String, u64>,
    pub unknown_bp_before: u64,
    pub unknown_bp_after: u64,
    /// Base pairs Unknown before, by their class after.
    pub unknown_resolved_to: BTreeMap<String, u64>,
    /// Base pairs Unknown before and labelled after, by what decided the label
    /// (`decided_by`; tandem-array satellites are cascade calls).
    pub unknown_resolved_by: BTreeMap<String, u64>,
    /// Families relabelled from `rep discover`.
    pub overridden_families: usize,
}

pub const RELABELLED_OUT: &str = "repeatmasker.relabelled.out";
pub const SUMMARY_TSV: &str = "annotation_summary.tsv";
pub const SUMMARY_JSON: &str = "annotation_summary.json";

/// The class of a RepeatMasker or `rep` label: `LINE/L2` → `LINE`.
fn class_of(label: &str) -> &str {
    label.split('/').next().unwrap_or(label)
}

fn is_unknown(label: &str) -> bool {
    class_of(label) == "Unknown"
}

/// Family name → (new label, decided_by), from `classification.tsv`.
fn read_labels(path: &Path) -> Result<HashMap<String, (String, String)>> {
    let mut labels = HashMap::new();
    let mut lines = BufReader::new(fs::File::open(path)?).lines();
    let header = lines.next().transpose()?.unwrap_or_default();
    let col = |name: &str| {
        header.split('\t').position(|c| c == name).ok_or_else(|| {
            Error::new(ErrorKind::GenericCli(format!(
                "{} has no {name} column",
                path.display()
            )))
        })
    };
    let (id, label) = (col("family_id")?, col("label")?);
    let decided = header.split('\t').position(|c| c == "decided_by");
    for line in lines {
        let line = line?;
        let f = line.split('\t').collect::<Vec<_>>();
        let (Some(name), Some(new)) = (f.get(id), f.get(label)) else {
            continue;
        };
        if new.is_empty() {
            continue;
        }
        let name = name.split('#').next().unwrap_or(name);
        let by = decided.and_then(|i| f.get(i)).copied().unwrap_or("cascade");
        labels.insert(name.to_owned(), ((*new).to_owned(), by.to_owned()));
    }
    Ok(labels)
}

/// Family name → label from `rep discover`'s `relabels.tsv`, for one sample.
fn read_overrides(path: &Path, sample: Option<&str>) -> Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    for line in BufReader::new(fs::File::open(path)?).lines().skip(1) {
        let line = line?;
        let f = line.split('\t').collect::<Vec<_>>();
        if f.len() < 4 || sample.is_some_and(|s| s != f[0]) {
            continue;
        }
        let name = f[1].split('#').next().unwrap_or(f[1]);
        out.push((name.to_owned(), f[3].to_owned()));
    }
    Ok(out)
}

/// Byte span of the `n`th whitespace-separated field of a line.
fn field_span(line: &str, n: usize) -> Option<(usize, usize)> {
    let mut seen = 0;
    let mut start = None;
    for (i, c) in line.char_indices() {
        match (c.is_whitespace(), start) {
            (false, None) => start = Some(i),
            (true, Some(s)) => {
                if seen == n {
                    return Some((s, i));
                }
                seen += 1;
                start = None;
            }
            _ => {}
        }
    }
    start.filter(|_| seen == n).map(|s| (s, line.len()))
}

/// Replace the class/family column, keeping the columns after it aligned
/// where the new label fits in the old one's padding.
fn rewrite(line: &str, span: (usize, usize), label: &str) -> String {
    let (s, e) = span;
    let pad_end = line[e..]
        .find(|c: char| !c.is_whitespace())
        .map_or(line.len(), |p| e + p);
    let width = pad_end - s;
    let cell = if label.len() < width {
        format!("{label:<width$}")
    } else {
        format!("{label} ")
    };
    format!("{}{}{}", &line[..s], cell, &line[pad_end..])
}

/// Total length of the union of half-open intervals.
fn union_bp(intervals: &mut [(u64, u64)]) -> u64 {
    intervals.sort_unstable();
    let (mut total, mut cur) = (0, None::<(u64, u64)>);
    for &(s, e) in intervals.iter() {
        cur = match cur {
            Some((cs, ce)) if s <= ce => Some((cs, ce.max(e))),
            Some((cs, ce)) => {
                total += ce - cs;
                Some((s, e))
            }
            None => Some((s, e)),
        };
    }
    total + cur.map_or(0, |(s, e)| e - s)
}

/// Base pairs per key, intervals keyed by (key, contig).
fn bp_by<K: Ord + Clone>(intervals: HashMap<(K, String), Vec<(u64, u64)>>) -> BTreeMap<K, u64> {
    let mut out = BTreeMap::new();
    for ((k, _), mut v) in intervals {
        *out.entry(k).or_default() += union_bp(&mut v);
    }
    out
}

fn read_sizes(path: &Path) -> Result<u64> {
    let mut total = 0;
    for line in BufReader::new(fs::File::open(path)?).lines() {
        let line = line?;
        if let Some(len) = line
            .split('\t')
            .nth(1)
            .and_then(|l| l.trim().parse::<u64>().ok())
        {
            total += len;
        }
    }
    Ok(total)
}

pub fn relabel(c: &RelabelConfig) -> Result<RelabelSummary> {
    fs::create_dir_all(&c.output_dir)?;
    let mut labels = read_labels(&c.classify_dir.join("classification.tsv"))?;
    let mut overridden = 0;
    if let Some(path) = &c.overrides {
        for (family, label) in read_overrides(path, c.sample.as_deref())? {
            labels.insert(family, (label, "discover".into()));
            overridden += 1;
        }
    }
    let mut out = BufWriter::new(fs::File::create(c.output_dir.join(RELABELLED_OUT))?);
    let mut s = RelabelSummary {
        genome_bp: c.genome_sizes.as_deref().map(read_sizes).transpose()?,
        overridden_families: overridden,
        ..Default::default()
    };
    type Intervals<K> = HashMap<(K, String), Vec<(u64, u64)>>;
    let (mut before, mut after): (Intervals<String>, Intervals<String>) = Default::default();
    let (mut all, mut unknown_before, mut unknown_after): (
        Intervals<()>,
        Intervals<()>,
        Intervals<()>,
    ) = Default::default();
    let (mut resolved_to, mut resolved_by): (Intervals<String>, Intervals<String>) =
        Default::default();
    for line in BufReader::new(fs::File::open(&c.rm_out)?).lines() {
        let line = line?;
        let f = line.split_whitespace().collect::<Vec<_>>();
        let hit = f.len() >= 15 && f[0].parse::<u64>().is_ok();
        let (Some(span), true) = (field_span(&line, 10), hit) else {
            writeln!(out, "{line}")?;
            continue;
        };
        let (Ok(start), Ok(end)) = (f[5].parse::<u64>(), f[6].parse::<u64>()) else {
            writeln!(out, "{line}")?;
            continue;
        };
        let iv = (start.saturating_sub(1), end);
        let contig = f[4].to_owned();
        let old = f[10];
        s.hits += 1;
        let (new, by) = match labels.get(f[9]) {
            Some((l, by)) if l != old => {
                s.relabelled_hits += 1;
                writeln!(out, "{}", rewrite(&line, span, l))?;
                (l.as_str(), by.as_str())
            }
            _ => {
                writeln!(out, "{line}")?;
                (old, "")
            }
        };
        let push = |m: &mut Intervals<String>, k: &str| {
            m.entry((k.to_owned(), contig.clone()))
                .or_default()
                .push(iv)
        };
        push(&mut before, class_of(old));
        push(&mut after, class_of(new));
        all.entry(((), contig.clone())).or_default().push(iv);
        if is_unknown(old) {
            unknown_before
                .entry(((), contig.clone()))
                .or_default()
                .push(iv);
            if !is_unknown(new) {
                push(&mut resolved_to, class_of(new));
                push(&mut resolved_by, if by.is_empty() { "cascade" } else { by });
            }
        }
        if is_unknown(new) {
            unknown_after
                .entry(((), contig.clone()))
                .or_default()
                .push(iv);
        }
    }
    out.flush()?;
    let total = |m: Intervals<()>| bp_by(m).into_values().sum::<u64>();
    s.masked_bp = total(all);
    s.unknown_bp_before = total(unknown_before);
    s.unknown_bp_after = total(unknown_after);
    s.bp_before = bp_by(before);
    s.bp_after = bp_by(after);
    s.unknown_resolved_to = bp_by(resolved_to);
    s.unknown_resolved_by = bp_by(resolved_by);

    let mut tsv = String::from("class\tbp_before\tbp_after\tpct_before\tpct_after\n");
    let classes = s
        .bp_before
        .keys()
        .chain(s.bp_after.keys())
        .collect::<std::collections::BTreeSet<_>>();
    let pct = |bp: u64| {
        s.genome_bp.map_or(String::new(), |g| {
            format!("{:.3}", 100.0 * bp as f64 / g.max(1) as f64)
        })
    };
    for class in classes {
        let (b, a) = (
            s.bp_before.get(class).copied().unwrap_or(0),
            s.bp_after.get(class).copied().unwrap_or(0),
        );
        tsv.push_str(&format!("{class}\t{b}\t{a}\t{}\t{}\n", pct(b), pct(a)));
    }
    fs::write(c.output_dir.join(SUMMARY_TSV), tsv)?;
    fs::write(
        c.output_dir.join(SUMMARY_JSON),
        serde_json::to_vec_pretty(&s)
            .map_err(|e| Error::new(ErrorKind::GenericCli(e.to_string())))?,
    )?;
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    const OUT: &str = "   SW   perc perc perc  query      position in query           matching       repeat              position in repeat
score   div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID

  1000   10.0  0.0  0.0  chr1            1      100    (900) + fam-1          Unknown                  1    100    (0)      1
  1000   10.0  0.0  0.0  chr1           51      150    (850) C fam-2          Unknown                (0)    100      1      2
   500   12.0  0.0  0.0  chr1          301      400    (600) + fam-3          LINE/L2                  1    100    (0)      3
    20    0.0  0.0  0.0  chr1          501      530    (470) + (TTAGG)n       Simple_repeat            1     30    (0)      4
";

    #[test]
    fn relabels_out_lines_and_counts_union_bp() {
        let dir = tempfile::tempdir().unwrap();
        let classify = dir.path().join("classify");
        fs::create_dir_all(&classify).unwrap();
        fs::write(
            classify.join("classification.tsv"),
            "family_id\toriginal_label\tlabel\tdecided_by\n\
             fam-1#Unknown\tUnknown\tDNA/hAT\trelative\n\
             fam-2#Unknown\tUnknown\tUnknown\tcascade\n\
             fam-3#LINE/L2\tLINE/L2\tLINE/RTE-BovB\tmodel\n",
        )
        .unwrap();
        let rm_out = dir.path().join("x.out");
        fs::write(&rm_out, OUT).unwrap();
        let sizes = dir.path().join("x.sizes");
        fs::write(&sizes, "chr1\t1000\n").unwrap();
        let s = relabel(&RelabelConfig {
            classify_dir: classify,
            rm_out,
            genome_sizes: Some(sizes),
            overrides: None,
            sample: None,
            output_dir: dir.path().join("out"),
        })
        .unwrap();
        assert_eq!((s.hits, s.relabelled_hits), (4, 2));
        assert_eq!(s.masked_bp, 150 + 100 + 30);
        assert_eq!(s.bp_before["Unknown"], 150);
        assert_eq!((s.bp_after["DNA"], s.bp_after["Unknown"]), (100, 100));
        assert_eq!((s.unknown_bp_before, s.unknown_bp_after), (150, 100));
        assert_eq!(s.unknown_resolved_to["DNA"], 100);
        assert_eq!(s.unknown_resolved_by["relative"], 100);
        assert_eq!(s.bp_after["LINE"], 100);

        let text = fs::read_to_string(dir.path().join("out").join(RELABELLED_OUT)).unwrap();
        let lines = text.lines().collect::<Vec<_>>();
        assert_eq!(lines.len(), OUT.lines().count());
        assert_eq!(lines[0], OUT.lines().next().unwrap());
        // same width, columns after the label stay aligned
        assert_eq!(lines[3].len(), OUT.lines().nth(3).unwrap().len());
        assert_eq!(lines[3].split_whitespace().nth(10), Some("DNA/hAT"));
        assert_eq!(lines[3].split_whitespace().nth(11), Some("1"));
        assert_eq!(lines[5].split_whitespace().nth(10), Some("LINE/RTE-BovB"));
        assert_eq!(lines[6], OUT.lines().nth(6).unwrap());
        let tsv = fs::read_to_string(dir.path().join("out").join(SUMMARY_TSV)).unwrap();
        assert!(tsv.contains("DNA\t0\t100\t0.000\t10.000\n"));
    }

    #[test]
    fn discover_relabels_override_classify_labels_for_their_sample() {
        let dir = tempfile::tempdir().unwrap();
        let classify = dir.path().join("classify");
        fs::create_dir_all(&classify).unwrap();
        fs::write(
            classify.join("classification.tsv"),
            "family_id\tlabel\tdecided_by\nfam-1#Unknown\tSINE/tRNA\trelative\nfam-2#Unknown\tUnknown\tcascade\n",
        )
        .unwrap();
        let rm_out = dir.path().join("x.out");
        fs::write(&rm_out, OUT).unwrap();
        let overrides = dir.path().join("relabels.tsv");
        fs::write(
            &overrides,
            "sample\tfamily_id\tlabel\tnew_label\tsource\tpurity\tlabelled_neighbours\n\
             sp/acc\tfam-1#Unknown\tSINE/tRNA\tDNA/hAT\tgroup\t1.00\t5\n\
             other/acc\tfam-2#Unknown\tUnknown\tLTR/Gypsy\tgroup\t1.00\t5\n",
        )
        .unwrap();
        let s = relabel(&RelabelConfig {
            classify_dir: classify,
            rm_out,
            genome_sizes: None,
            overrides: Some(overrides),
            sample: Some("sp/acc".into()),
            output_dir: dir.path().join("out"),
        })
        .unwrap();
        assert_eq!(s.overridden_families, 1);
        assert_eq!(s.unknown_resolved_by["discover"], 100);
        assert_eq!(s.unknown_resolved_to["DNA"], 100);
        assert_eq!(s.bp_after["Unknown"], 100);
        assert!(!s.bp_after.contains_key("SINE"));
    }

    #[test]
    fn a_label_longer_than_its_column_keeps_fields_apart() {
        let line = "  1 0 0 0 c 1 2 (0) + f Unknown 1 2 (0) 1";
        let span = field_span(line, 10).unwrap();
        let new = rewrite(line, span, "DNA/CMC-EnSpm");
        assert_eq!(
            new.split_whitespace().collect::<Vec<_>>()[9..12],
            ["f", "DNA/CMC-EnSpm", "1"]
        );
    }
}
