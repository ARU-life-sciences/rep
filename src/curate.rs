//! `rep curate`: refined consensi from the genomic copies `rep recover-copies` found.
//!
//! For each family the copies, aligned with their flanks, are anchored on the
//! columns where they matched the query (their BLAST cores). Boundaries only
//! move outward from that anchor, while the copies keep agreeing, and stop
//! where agreement falls below a floor or drops sharply: a block of copies
//! ending together, such as solo LTRs beside the LTRs of full elements. The
//! consensus is the majority base of the copies spanning each column. Columns
//! where most copies have a gap (insertions in a few copies) are dropped, never
//! written as N. A family keeps its query when it has too few copies, or when
//! the refined consensus would be much shorter than the query.
//!
//! The output has the layout of a `rep recover-copies` directory, with copy
//! boundaries moved to the refined ends, so `rep classify --copies` calls
//! full-length copies and TSDs against the refined consensus.

use std::{
    collections::HashMap,
    fs,
    path::{Path, PathBuf},
    sync::Mutex,
};

use bio::io::fasta;
use serde::Serialize;

use crate::{
    classify::structure::{self, StructureParams},
    copy_recovery::family_dir_name,
    Error, ErrorKind, Result,
};

pub const REFINED: &str = "refined.fa";
pub const TABLE: &str = "curation.tsv";
pub const SUMMARY: &str = "curation.json";

#[derive(Clone, Debug, Serialize)]
pub struct CurateParams {
    /// Consensus columns averaged when testing copy agreement.
    pub window: usize,
    /// Extension stops where the fraction of spanning copies sharing the
    /// majority base falls below this.
    pub agreement_floor: f64,
    /// ...or falls by this much from the preceding window.
    pub agreement_drop: f64,
    /// Fewest copies sharing a base for a column to count as agreeing.
    pub min_agreeing: usize,
    /// Fewest aligned copies needed to refine a family.
    pub min_copies: usize,
    /// Fraction of copies whose core must cover a column for it to be anchor.
    pub anchor_fraction: f64,
    /// A refined consensus shorter than this fraction of the query is rejected.
    pub min_length_fraction: f64,
    /// Shortest inverted terminal repeat the extended ends snap to, and its
    /// identity. Stricter than `rep classify`: a snapped end would otherwise
    /// manufacture the TIR evidence classify then reads.
    pub snap_min_tir: usize,
    pub snap_min_identity: f64,
    /// Extended ends are trimmed back by up to `tsd_max_trim` bases where
    /// that leaves a target-site duplication (up to `tsd_max_length` bp)
    /// outside at least `tsd_min_copies` copies and `tsd_min_support` of them.
    pub tsd_max_trim: usize,
    pub tsd_max_length: usize,
    pub tsd_min_support: f64,
    pub tsd_min_copies: usize,
    /// Terminal-repeat search windows and identity (shared with `rep classify`).
    pub structure: StructureParams,
}

impl Default for CurateParams {
    fn default() -> Self {
        Self {
            window: 25,
            agreement_floor: 0.5,
            agreement_drop: 0.25,
            min_agreeing: 3,
            min_copies: 3,
            anchor_fraction: 0.2,
            min_length_fraction: 0.8,
            snap_min_tir: 25,
            snap_min_identity: 0.9,
            tsd_max_trim: 10,
            tsd_max_length: 10,
            tsd_min_support: 0.7,
            tsd_min_copies: 5,
            structure: StructureParams::default(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct CurateConfig {
    /// `rep recover-copies` output directory.
    pub copies_dir: PathBuf,
    pub out: PathBuf,
    pub threads: usize,
    pub params: CurateParams,
}

impl CurateConfig {
    pub fn new(copies_dir: PathBuf, out: PathBuf) -> Self {
        Self {
            copies_dir,
            out,
            threads: 1,
            params: CurateParams::default(),
        }
    }
}

/// What happened to one family.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct FamilyCuration {
    pub family_id: String,
    pub query_length: usize,
    pub refined_length: usize,
    /// Aligned copies used.
    pub copies: usize,
    /// Consensus bases added beyond the query's aligned region, 5' and 3'.
    pub added_5: usize,
    pub added_3: usize,
    /// Copies spanning at least 80% of the refined consensus.
    pub full_length_copies: usize,
    /// `refined`, or why the query was kept.
    pub status: String,
}

#[derive(Debug, Default, Serialize)]
pub struct CurateSummary {
    pub families: usize,
    pub refined: usize,
    pub kept: HashMap<String, usize>,
    pub bases_added: usize,
}

/// One aligned copy.
#[derive(Clone, Debug)]
struct Row {
    locus_id: String,
    seq: Vec<u8>,
    /// First and last non-gap columns.
    first: usize,
    last: usize,
    /// Columns of the BLAST core (the copy minus its flanks).
    core: (usize, usize),
}

#[derive(Clone, Copy, Debug, Default)]
struct Column {
    /// Copies whose sequence spans the column.
    span: usize,
    gaps: usize,
    bases: usize,
    major: u8,
    major_count: usize,
}

fn column(rows: &[Row], c: usize, core_only: bool) -> Column {
    let mut counts = [0usize; 4];
    let mut col = Column::default();
    for r in rows {
        if c < r.first || c > r.last || (core_only && (c < r.core.0 || c > r.core.1)) {
            continue;
        }
        col.span += 1;
        match r.seq[c] {
            b'-' | b'.' => col.gaps += 1,
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            _ => {}
        }
    }
    col.bases = counts.iter().sum();
    let (i, &n) = counts
        .iter()
        .enumerate()
        .max_by_key(|&(i, n)| (*n, std::cmp::Reverse(i)))
        .expect("four bases");
    col.major = b"ACGT"[i];
    col.major_count = n;
    col
}

/// The refined consensus of one family, or why there is none.
#[derive(Debug, PartialEq)]
struct Refined {
    seq: Vec<u8>,
    /// Alignment column of each consensus base.
    columns: Vec<usize>,
    added_5: usize,
    added_3: usize,
}

fn refine(rows: &[Row], query_len: usize, p: &CurateParams) -> std::result::Result<Refined, &'static str> {
    if rows.len() < p.min_copies {
        return Err("too_few_copies");
    }
    let width = rows[0].seq.len();
    // copies whose core covers each column
    let mut steps = vec![0i64; width + 1];
    for r in rows {
        steps[r.core.0] += 1;
        steps[r.core.1 + 1] -= 1;
    }
    let core_n = steps
        .iter()
        .scan(0i64, |n, s| {
            *n += s;
            Some(*n as usize)
        })
        .collect::<Vec<_>>();
    let need = 2.max((p.anchor_fraction * rows.len() as f64).ceil() as usize);
    let (Some(a0), Some(a1)) = (
        (0..width).find(|&c| core_n[c] >= need),
        (0..width).rev().find(|&c| core_n[c] >= need),
    ) else {
        return Err("no_anchor");
    };
    let all = (0..width).map(|c| column(rows, c, false)).collect::<Vec<_>>();
    // agreement along consensus columns (minority-insertion columns skipped)
    let agreement = |cols: &mut dyn Iterator<Item = usize>| {
        let mut kept = Vec::new();
        let mut agree = Vec::new();
        for c in cols {
            let s = all[c];
            if s.span == 0 {
                break;
            }
            if s.gaps > s.bases {
                continue;
            }
            kept.push(c);
            agree.push(if s.major_count >= p.min_agreeing {
                s.major_count as f64 / s.span as f64
            } else {
                0.0
            });
        }
        (kept, agree)
    };
    let w = p.window.max(1);
    // walk outward from the anchor edge; `inner` ends at the edge, `outer` leaves it
    let extend = |inner: &mut dyn Iterator<Item = usize>, outer: &mut dyn Iterator<Item = usize>| {
        let (_, inside) = agreement(inner);
        let (cols, outside) = agreement(outer);
        let inside = &inside[inside.len().saturating_sub(w)..];
        let a = [inside, &outside[..]].concat();
        let k0 = inside.len();
        let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        for k in k0..a.len() {
            let after = &a[k..(k + w).min(a.len())];
            let before = &a[k.saturating_sub(w)..k];
            let before = if before.is_empty() { 1.0 } else { mean(before) };
            let after = mean(after);
            if after < p.agreement_floor || before - after >= p.agreement_drop {
                // the windows straddle the end: place it at the best two-level split
                let (lo, hi) = (k.saturating_sub(w).max(1), (k + 2 * w).min(a.len()));
                let sse = |v: &[f64]| {
                    let m = mean(v);
                    v.iter().map(|x| (x - m).powi(2)).sum::<f64>()
                };
                let end = (k.max(k0)..=(k + w).min(hi))
                    .filter(|&j| j > lo && j < hi)
                    .min_by(|&x, &y| {
                        (sse(&a[lo..x]) + sse(&a[x..hi])).total_cmp(&(sse(&a[lo..y]) + sse(&a[y..hi])))
                    })
                    .unwrap_or(k);
                return cols[..end.max(k0) - k0].to_vec();
            }
        }
        cols
    };
    let reach = 3 * w;
    let right = extend(
        &mut (a1.saturating_sub(reach).max(a0)..=a1),
        &mut (a1 + 1..width),
    );
    let left = extend(
        &mut (a0..=(a0 + reach).min(a1)).rev(),
        &mut (0..a0).rev(),
    );
    let lo = left.last().copied().unwrap_or(a0);
    let hi = right.last().copied().unwrap_or(a1);
    let mut seq = Vec::new();
    let mut columns = Vec::new();
    for c in lo..=hi {
        // inside the anchor, copies whose core does not reach a column are flank there
        let s = if c >= a0 && c <= a1 && core_n[c] >= 2 {
            column(rows, c, true)
        } else {
            all[c]
        };
        if s.span == 0 || s.gaps > s.bases || s.bases == 0 {
            continue;
        }
        seq.push(s.major);
        columns.push(c);
    }
    let mut added_5 = columns.iter().take_while(|&&c| c < a0).count();
    let mut added_3 = columns.iter().rev().take_while(|&&c| c > a1).count();
    let (trim_5, trim_3) = snap_to_terminal_repeat(&seq, added_5, added_3, p);
    seq.truncate(seq.len() - trim_3);
    seq.drain(..trim_5);
    columns.truncate(columns.len() - trim_3);
    columns.drain(..trim_5);
    added_5 -= trim_5;
    added_3 -= trim_3;
    let (trim_5, trim_3) = trim_to_tsd(rows, &columns, added_5, added_3, p);
    seq.truncate(seq.len() - trim_3);
    seq.drain(..trim_5);
    columns.truncate(columns.len() - trim_3);
    columns.drain(..trim_5);
    added_5 -= trim_5;
    added_3 -= trim_3;
    if (seq.len() as f64) < p.min_length_fraction * query_len as f64 {
        return Err("shorter_than_query");
    }
    Ok(Refined {
        seq,
        columns,
        added_5,
        added_3,
    })
}

/// Bases to trim from the extended ends so that the consensus ends with a
/// terminal repeat (direct, LTR-like, tried first; then inverted, TIR-like)
/// whose outer edges lie in the added sequence. Copies of an element nested in
/// or next to another repeat can agree past its true ends; its own terminal
/// repeats mark them.
fn snap_to_terminal_repeat(
    seq: &[u8],
    added_5: usize,
    added_3: usize,
    p: &CurateParams,
) -> (usize, usize) {
    if added_5 + added_3 == 0 {
        return (0, 0);
    }
    let slack = added_5.max(added_3);
    let s = &p.structure;
    for (inverted, window, min_length, identity) in [
        (false, s.ltr_window, s.ltr_min_length, s.terminal_min_identity),
        (true, s.tir_window, p.snap_min_tir, p.snap_min_identity),
    ] {
        let window = window + slack;
        if let Some(r) = structure::terminal_repeat(
            seq,
            inverted,
            window,
            min_length,
            identity,
            slack,
        ) {
            if r.five_prime_offset <= added_5 && r.three_prime_offset <= added_3 {
                return (r.five_prime_offset, r.three_prime_offset);
            }
        }
    }
    (0, 0)
}

/// Bases to trim from the extended ends so that the copies' flanks carry a
/// target-site duplication just outside them. A TSD that is the same in
/// every copy (TA for Tc1/mariner) agrees as well as the element does, so the
/// extension takes it in; copies' flanks show where the element really ends.
fn trim_to_tsd(
    rows: &[Row],
    columns: &[usize],
    added_5: usize,
    added_3: usize,
    p: &CurateParams,
) -> (usize, usize) {
    // bases beyond a duplication, which must differ between copies: within
    // the element, both sides agree everywhere and chance 2-mer matches pass
    const OUTSIDE: usize = 4;
    let n = columns.len();
    let k = p.tsd_max_length;
    if added_5 + added_3 == 0 || n < 2 * p.tsd_max_trim + 2 || k < 2 {
        return (0, 0);
    }
    let ungapped = |s: &[u8]| s.iter().copied().filter(|&b| b != b'-' && b != b'.').collect::<Vec<_>>();
    // flanks of the copies spanning both ends when trimmed by (t5, t3)
    let flanks = |t5: usize, t3: usize| {
        let (cs, ce) = (columns[t5], columns[n - 1 - t3]);
        rows.iter()
            .filter(|r| r.first < cs && r.last > ce)
            .map(|r| (ungapped(&r.seq[..cs]), ungapped(&r.seq[ce + 1..])))
            .filter(|(l, r)| l.len() >= k + OUTSIDE && r.len() >= k + OUTSIDE)
            .collect::<Vec<_>>()
    };
    let has = |(l, r): &(Vec<u8>, Vec<u8>), m: usize| l[l.len() - m..] == r[..m];
    let fraction = |f: &[(Vec<u8>, Vec<u8>)], m: usize| {
        if f.is_empty() {
            0.0
        } else {
            f.iter().filter(|x| has(x, m)).count() as f64 / f.len() as f64
        }
    };
    // share of the majority base just outside the duplication, on the more
    // conserved side
    let outside_agreement = |f: &[&(Vec<u8>, Vec<u8>)], m: usize| {
        [true, false]
            .into_iter()
            .map(|left| {
                let mut sum = 0.0;
                for j in 0..OUTSIDE {
                    let mut counts = HashMap::new();
                    for (l, r) in f {
                        let b = if left { l[l.len() - m - 1 - j] } else { r[m + j] };
                        *counts.entry(b).or_insert(0usize) += 1;
                    }
                    sum += *counts.values().max().unwrap_or(&0) as f64 / f.len() as f64;
                }
                sum / OUTSIDE as f64
            })
            .fold(0.0, f64::max)
    };
    let base = flanks(0, 0);
    let mut best: Option<(f64, usize, usize)> = None;
    for t5 in 0..=added_5.min(p.tsd_max_trim) {
        for t3 in 0..=added_3.min(p.tsd_max_trim) {
            if t5 + t3 == 0 {
                continue;
            }
            let f = flanks(t5, t3);
            let Some(m) = (2..=k).rev().find(|&m| {
                let with = f.iter().filter(|x| has(x, m)).count();
                with >= p.tsd_min_copies && with as f64 >= p.tsd_min_support * f.len() as f64
            }) else {
                continue;
            };
            let frac = fraction(&f, m);
            let with = f.iter().filter(|x| has(x, m)).collect::<Vec<_>>();
            if frac < 2.0 * fraction(&base, m) || outside_agreement(&with, m) >= 0.75 {
                continue;
            }
            if best.is_none_or(|(b, a5, a3)| frac > b || (frac == b && t5 + t3 < a5 + a3)) {
                best = Some((frac, t5, t3));
            }
        }
    }
    best.map_or((0, 0), |(_, t5, t3)| (t5, t3))
}

/// Where a copy matches the refined consensus: its first and last consensus
/// positions (where at least half of the next `window` positions match), and
/// its identity over that span.
fn copy_span(row: &Row, r: &Refined, window: usize) -> Option<(usize, usize, f64)> {
    let m = r
        .columns
        .iter()
        .zip(&r.seq)
        .map(|(&c, &b)| row.seq[c] == b)
        .collect::<Vec<_>>();
    let w = window.min(m.len()).max(1);
    let good = |k: usize, back: bool| {
        let (a, z) = if back {
            ((k + 1).saturating_sub(w), k + 1)
        } else {
            (k, (k + w).min(m.len()))
        };
        m[k] && 2 * m[a..z].iter().filter(|&&x| x).count() >= z - a
    };
    let start = (0..m.len()).find(|&k| good(k, false))?;
    let end = (start..m.len()).rev().find(|&k| good(k, true))?;
    let same = m[start..=end].iter().filter(|&&x| x).count();
    Some((start, end, 100.0 * same as f64 / (end - start + 1) as f64))
}

/// Ungapped bases of `row` before alignment column `c`.
fn bases_before(row: &Row, c: usize) -> usize {
    row.seq[..c].iter().filter(|&&b| b != b'-' && b != b'.').count()
}

fn bases_after(row: &Row, c: usize) -> usize {
    row.seq[c + 1..]
        .iter()
        .filter(|&&b| b != b'-' && b != b'.')
        .count()
}

/// New position of one copy against the refined consensus.
#[derive(Debug)]
struct CopyFit {
    /// Flanks in consensus orientation.
    left: usize,
    right: usize,
    coverage: f64,
    identity: f64,
}

fn read_fasta(path: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    fasta::Reader::from_file(path)
        .map_err(|e| Error::new(ErrorKind::GenericCli(format!("{}: {e}", path.display()))))?
        .records()
        .map(|r| {
            let r = r?;
            Ok((r.id().to_owned(), r.seq().to_ascii_uppercase()))
        })
        .collect()
}

/// Oriented flanks of the selected copies, from `copies.tsv`.
fn read_flanks(path: &Path) -> Result<HashMap<String, (usize, usize)>> {
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    let headers = reader.headers()?.clone();
    let col = |name: &str| headers.iter().position(|h| h == name);
    let (Some(id), Some(left), Some(right), Some(rc)) = (
        col("locus_id"),
        col("left_flank_actual"),
        col("right_flank_actual"),
        col("reverse_complemented"),
    ) else {
        return Ok(HashMap::new());
    };
    let mut flanks = HashMap::new();
    for row in reader.records() {
        let row = row?;
        let (l, r): (usize, usize) = (row[left].parse()?, row[right].parse()?);
        // flanks are in genomic orientation; minus-strand copies were reverse-complemented
        flanks.insert(
            row[id].to_owned(),
            if &row[rc] == "true" { (r, l) } else { (l, r) },
        );
    }
    Ok(flanks)
}

fn rows_of(dir: &Path) -> Result<Vec<Row>> {
    let aligned = dir.join("copies.aligned.fa");
    if !aligned.is_file() {
        return Ok(Vec::new());
    }
    let flanks = read_flanks(&dir.join("copies.tsv"))?;
    let mut rows = Vec::new();
    for (id, seq) in read_fasta(&aligned)? {
        let locus = id.trim_start_matches("locus_id=").to_owned();
        let Some(&(left, right)) = flanks.get(&locus) else {
            continue;
        };
        let cols = seq
            .iter()
            .enumerate()
            .filter(|(_, &b)| b != b'-' && b != b'.')
            .map(|(i, _)| i)
            .collect::<Vec<_>>();
        if cols.len() <= left + right {
            continue;
        }
        rows.push(Row {
            locus_id: locus,
            first: cols[0],
            last: cols[cols.len() - 1],
            core: (cols[left], cols[cols.len() - 1 - right]),
            seq,
        });
    }
    if rows.windows(2).any(|w| w[0].seq.len() != w[1].seq.len()) {
        return Err(Error::new(ErrorKind::GenericCli(format!(
            "{}: aligned copies differ in length",
            aligned.display()
        ))));
    }
    Ok(rows)
}

/// Rewrite a recover-copies TSV, applying `edit` to rows of fitted copies.
fn rewrite_tsv(
    from: &Path,
    to: &Path,
    mut edit: impl FnMut(&HashMap<&str, usize>, &mut Vec<String>),
) -> Result<()> {
    if !from.is_file() {
        return Ok(());
    }
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(from)?;
    let headers = reader.headers()?.clone();
    let index = headers
        .iter()
        .enumerate()
        .map(|(i, h)| (h, i))
        .collect::<HashMap<_, _>>();
    let mut out = headers.iter().collect::<Vec<_>>().join("\t") + "\n";
    for row in reader.records() {
        let mut fields = row?.iter().map(str::to_owned).collect::<Vec<_>>();
        edit(&index, &mut fields);
        out += &(fields.join("\t") + "\n");
    }
    fs::write(to, out)?;
    Ok(())
}

fn curate_family(dir: &Path, out: &Path, p: &CurateParams) -> Result<(FamilyCuration, Vec<u8>)> {
    let query = read_fasta(&dir.join("query.fa"))?
        .into_iter()
        .next()
        .ok_or_else(|| Error::new(ErrorKind::GenericCli(format!("{}: no query", dir.display()))))?;
    let (family_id, query_seq) = query;
    let rows = rows_of(dir)?;
    let fam_out = out.join("families").join(family_dir_name(&family_id));
    fs::create_dir_all(&fam_out)?;
    let mut record = FamilyCuration {
        family_id: family_id.clone(),
        query_length: query_seq.len(),
        refined_length: query_seq.len(),
        copies: rows.len(),
        added_5: 0,
        added_3: 0,
        full_length_copies: 0,
        status: String::new(),
    };
    let extended = dir.join("copies.extended.fa");
    if extended.is_file() {
        fs::copy(&extended, fam_out.join("copies.extended.fa"))?;
    }
    let refined = match refine(&rows, query_seq.len(), p) {
        Ok(r) => r,
        Err(why) => {
            record.status = why.to_owned();
            for f in ["copies.tsv", "loci.tsv"] {
                if dir.join(f).is_file() {
                    fs::copy(dir.join(f), fam_out.join(f))?;
                }
            }
            fs::write(
                fam_out.join("query.fa"),
                format!(">{family_id}\n{}\n", String::from_utf8_lossy(&query_seq)),
            )?;
            return Ok((record, query_seq));
        }
    };
    let mut fits = HashMap::new();
    for row in &rows {
        if let Some((s, e, identity)) = copy_span(row, &refined, p.window) {
            let coverage = (e - s + 1) as f64 / refined.seq.len() as f64;
            fits.insert(
                row.locus_id.clone(),
                CopyFit {
                    left: bases_before(row, refined.columns[s]),
                    right: bases_after(row, refined.columns[e]),
                    coverage,
                    identity,
                },
            );
        }
    }
    record.status = "refined".into();
    record.refined_length = refined.seq.len();
    record.added_5 = refined.added_5;
    record.added_3 = refined.added_3;
    record.full_length_copies = fits.values().filter(|f| f.coverage >= 0.8).count();

    // copies.tsv: flanks and genomic ends of fitted copies at the refined boundaries
    rewrite_tsv(&dir.join("copies.tsv"), &fam_out.join("copies.tsv"), |ix, f| {
        let Some(fit) = fits.get(&f[ix["locus_id"]]) else {
            return;
        };
        let (left, right) = if f[ix["reverse_complemented"]] == "true" {
            (fit.right, fit.left)
        } else {
            (fit.left, fit.right)
        };
        let es: u64 = f[ix["actual_extended_start"]].parse().unwrap_or(0);
        let ee: u64 = f[ix["actual_extended_end"]].parse().unwrap_or(0);
        f[ix["left_flank_actual"]] = left.to_string();
        f[ix["right_flank_actual"]] = right.to_string();
        f[ix["genomic_start"]] = (es + left as u64).to_string();
        f[ix["genomic_end"]] = ee.saturating_sub(right as u64).to_string();
    })?;
    // loci.tsv: coverage of the refined consensus; loci not aligned are rescaled
    let scale = query_seq.len() as f64 / refined.seq.len() as f64;
    rewrite_tsv(&dir.join("loci.tsv"), &fam_out.join("loci.tsv"), |ix, f| {
        let (Some(&cov), Some(&id)) = (ix.get("combined_query_coverage"), ix.get("percent_identity")) else {
            return;
        };
        if let Some(fit) = fits.get(&f[ix["locus_id"]]) {
            f[cov] = format!("{:.4}", fit.coverage);
            f[id] = format!("{:.3}", fit.identity);
        } else if let Ok(c) = f[cov].parse::<f64>() {
            f[cov] = format!("{:.4}", (c * scale).min(1.0));
        }
    })?;
    fs::write(
        fam_out.join("query.fa"),
        format!(">{family_id}\n{}\n", String::from_utf8_lossy(&refined.seq)),
    )?;
    Ok((record, refined.seq))
}

/// Refine every family of a `rep recover-copies` directory.
pub fn run_curate(c: &CurateConfig) -> Result<CurateSummary> {
    let families = c.copies_dir.join("families");
    let mut dirs = fs::read_dir(&families)
        .map_err(|e| Error::new(ErrorKind::GenericCli(format!("{}: {e}", families.display()))))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|d| d.join("query.fa").is_file())
        .collect::<Vec<_>>();
    dirs.sort();
    fs::create_dir_all(c.out.join("families"))?;
    let results = Mutex::new(Vec::new());
    let next = Mutex::new(0usize);
    std::thread::scope(|scope| {
        for _ in 0..c.threads.max(1) {
            scope.spawn(|| loop {
                let i = {
                    let mut n = next.lock().expect("lock");
                    let i = *n;
                    *n += 1;
                    i
                };
                let Some(dir) = dirs.get(i) else {
                    break;
                };
                let r = curate_family(dir, &c.out, &c.params);
                results.lock().expect("lock").push((i, r));
            });
        }
    });
    let mut results = results.into_inner().expect("lock");
    results.sort_by_key(|(i, _)| *i);

    let mut summary = CurateSummary::default();
    let mut fasta = String::new();
    let mut table = String::from(
        "family_id\tquery_length\trefined_length\tcopies\tadded_5\tadded_3\tfull_length_copies\tstatus\n",
    );
    for (_, r) in results {
        let (f, seq) = r?;
        summary.families += 1;
        if f.status == "refined" {
            summary.refined += 1;
            summary.bases_added += f.added_5 + f.added_3;
        } else {
            *summary.kept.entry(f.status.clone()).or_default() += 1;
        }
        fasta += &format!(">{}\n", f.family_id);
        for line in seq.chunks(60) {
            fasta += std::str::from_utf8(line).expect("ASCII");
            fasta.push('\n');
        }
        table += &format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            f.family_id,
            f.query_length,
            f.refined_length,
            f.copies,
            f.added_5,
            f.added_3,
            f.full_length_copies,
            f.status
        );
    }
    fs::write(c.out.join(REFINED), fasta)?;
    fs::write(c.out.join(TABLE), table)?;
    fs::write(
        c.out.join(SUMMARY),
        serde_json::to_string_pretty(&serde_json::json!({
            "copies_dir": c.copies_dir,
            "params": c.params,
            "summary": summary,
        }))
        .map_err(|e| Error::new(ErrorKind::GenericCli(e.to_string())))?,
    )?;
    Ok(summary)
}

#[cfg(test)]
mod tests;
