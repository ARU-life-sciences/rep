//! Copy evidence from an existing RepeatMasker annotation (`.out`).
//!
//! RepeatMasker already reports every genomic copy of every library family,
//! with divergence from the consensus and the covered consensus interval, and
//! links fragments of one insertion through its ID column. That gives copy
//! number, completeness and age without running copy recovery. With the
//! genome sequence, flanks of full-length copies are cut out by coordinate
//! for TSD and flank-conservation analysis.

use std::{
    collections::HashMap,
    fs,
    io::{BufRead, BufReader, Read},
    path::Path,
    process::{Command, Stdio},
};

use bio::{alphabets::dna, io::fasta};

use super::copies::{
    call_tsd, flank_similarity, helitron_ends, CopyEvidence, CopyParams, OrientedCopy,
};
use crate::{Error, ErrorKind, Result};

/// One line of a RepeatMasker `.out` file.
#[derive(Clone, Debug, PartialEq)]
pub struct RmHit {
    pub contig: String,
    /// 0-based half-open genomic interval.
    pub start: u64,
    pub end: u64,
    pub minus: bool,
    pub family: String,
    pub class: String,
    /// 1-based inclusive consensus interval.
    pub rep_start: u64,
    pub rep_end: u64,
    pub divergence: f64,
    pub score: u64,
    pub id: u64,
}

fn parse_line(line: &str) -> Option<RmHit> {
    let f = line.split_whitespace().collect::<Vec<_>>();
    if f.len() < 15 || f[0].parse::<u64>().is_err() {
        return None;
    }
    let minus = f[8] == "C";
    // consensus coordinates: `begin end (left)` on +, `(left) end begin` on C
    let (rep_start, rep_end) = if minus {
        (f[13].parse().ok()?, f[12].parse().ok()?)
    } else {
        (f[11].parse().ok()?, f[12].parse().ok()?)
    };
    Some(RmHit {
        contig: f[4].to_owned(),
        start: f[5].parse::<u64>().ok()?.saturating_sub(1),
        end: f[6].parse().ok()?,
        minus,
        family: f[9].to_owned(),
        class: f[10].to_owned(),
        rep_start,
        rep_end,
        divergence: f[1].parse().ok()?,
        score: f[0].parse().ok()?,
        id: f[14].parse().ok()?,
    })
}

/// Read a RepeatMasker `.out` file, grouped by family name.
pub fn read_rm_out(path: &Path) -> Result<HashMap<String, Vec<RmHit>>> {
    let mut families: HashMap<String, Vec<RmHit>> = HashMap::new();
    for line in BufReader::new(fs::File::open(path)?).lines() {
        if let Some(hit) = parse_line(&line?) {
            families.entry(hit.family.clone()).or_default().push(hit);
        }
    }
    Ok(families)
}

/// Load a genome FASTA (plain or gzip-compressed) into memory.
pub fn load_genome(path: &Path) -> Result<HashMap<String, Vec<u8>>> {
    let gz = path.extension().is_some_and(|e| e == "gz");
    let mut child = None;
    let reader: Box<dyn Read> = if gz {
        let mut c = Command::new("gzip")
            .arg("-dc")
            .arg(path)
            .stdout(Stdio::piped())
            .spawn()
            .map_err(|e| Error::new(ErrorKind::GenericCli(format!("could not run gzip: {e}"))))?;
        let out = c.stdout.take().expect("piped");
        child = Some(c);
        Box::new(out)
    } else {
        Box::new(fs::File::open(path)?)
    };
    let mut genome = HashMap::new();
    for record in fasta::Reader::new(reader).records() {
        let record = record?;
        genome.insert(record.id().to_owned(), record.seq().to_ascii_uppercase());
    }
    if let Some(mut c) = child {
        let status = c.wait()?;
        if !status.success() {
            return Err(Error::new(ErrorKind::GenericCli(format!(
                "gzip failed on {}",
                path.display()
            ))));
        }
    }
    Ok(genome)
}

/// One insertion: RepeatMasker fragments sharing an ID.
struct Insertion<'a> {
    hits: Vec<&'a RmHit>,
    coverage: f64,
    span: u64,
}

/// Flank length cut from the genome on each side, and core bases kept at each end.
/// Flanks leave room for the TSD null after extending copies to the consensus ends.
const FLANK: u64 = 70;
const CORE_END: u64 = 30;
/// Most full-length copies used for TSD and flank analysis per family.
const MAX_COPIES: usize = 100;

/// Copy evidence for one family from its RepeatMasker hits.
pub fn copy_evidence_from_rm(
    hits: &[RmHit],
    consensus_length: usize,
    genome: Option<&HashMap<String, Vec<u8>>>,
    params: &CopyParams,
) -> CopyEvidence {
    let len = consensus_length.max(1) as u64;
    let mut by_id: HashMap<u64, Vec<&RmHit>> = HashMap::new();
    for h in hits {
        by_id.entry(h.id).or_default().push(h);
    }
    let mut insertions = by_id
        .into_values()
        .map(|mut hits| {
            hits.sort_by_key(|h| h.rep_start);
            let mut covered = 0;
            let mut last = 0;
            for h in &hits {
                let (s, e) = (h.rep_start.saturating_sub(1).max(last), h.rep_end.min(len));
                if e > s {
                    covered += e - s;
                    last = e;
                }
            }
            let start = hits.iter().map(|h| h.start).min().unwrap_or(0);
            let end = hits.iter().map(|h| h.end).max().unwrap_or(0);
            Insertion {
                coverage: covered as f64 / len as f64,
                span: end - start,
                hits,
            }
        })
        .collect::<Vec<_>>();
    // deterministic order: best-scoring insertions first
    insertions.sort_by_key(|i| std::cmp::Reverse(i.hits.iter().map(|h| h.score).sum::<u64>()));

    let (mut weighted, mut total) = (0.0, 0.0);
    for h in hits {
        let l = (h.end - h.start) as f64;
        weighted += h.divergence * l;
        total += l;
    }
    // a genomic span far longer than the consensus is a tandem array, not a copy
    let is_full =
        |i: &Insertion| i.coverage >= params.full_length_coverage && i.span <= len * 3 / 2;
    let mut evidence = CopyEvidence {
        loci: insertions.len(),
        full_length: insertions.iter().filter(|i| is_full(i)).count(),
        mean_identity: if total > 0.0 {
            100.0 - weighted / total
        } else {
            0.0
        },
        tsd: None,
        flank_similarity: None,
        genome_bp: hits.iter().map(|h| h.end - h.start).sum(),
        helitron_ends: None,
    };

    if let Some(genome) = genome {
        let mut copies = Vec::new();
        for ins in insertions.iter().filter(|i| is_full(i)) {
            // the copy must reach both consensus ends for its boundaries to be the element's
            let reaches_start = ins.hits.iter().any(|h| h.rep_start <= 10);
            let reaches_end = ins.hits.iter().any(|h| h.rep_end + 10 >= len);
            if !(reaches_start && reaches_end) {
                continue;
            }
            let first = ins.hits[0];
            let minus = ins.hits.iter().filter(|h| h.minus).count() * 2 > ins.hits.len();
            let (start, end) = (
                ins.hits.iter().map(|h| h.start).min().unwrap_or(0),
                ins.hits.iter().map(|h| h.end).max().unwrap_or(0),
            );
            // extend to the consensus ends, so that all copies align at both
            let missing_start = ins.hits.iter().map(|h| h.rep_start).min().unwrap_or(1) - 1;
            let missing_end =
                len.saturating_sub(ins.hits.iter().map(|h| h.rep_end).max().unwrap_or(len));
            let (start, end) = if minus {
                (start.saturating_sub(missing_end), end + missing_start)
            } else {
                (start.saturating_sub(missing_start), end + missing_end)
            };
            if let Some(copy) = cut_copy(genome, &first.contig, start, end, minus) {
                copies.push(copy);
            }
            if copies.len() >= MAX_COPIES {
                break;
            }
        }
        let refs = copies.iter().collect::<Vec<_>>();
        evidence.tsd = call_tsd(&refs, params);
        evidence.helitron_ends = helitron_ends(&refs);
        evidence.flank_similarity = flank_similarity(&copies, params.flank_window);
    }
    evidence
}

/// Flanks plus the first and last bases of the copy, in consensus orientation.
fn cut_copy(
    genome: &HashMap<String, Vec<u8>>,
    contig: &str,
    start: u64,
    end: u64,
    minus: bool,
) -> Option<OrientedCopy> {
    let seq = genome.get(contig)?;
    if start < FLANK || end + FLANK > seq.len() as u64 || end <= start {
        return None;
    }
    let core = end - start;
    let (s, e) = (start as usize, end as usize);
    let (f, c) = (FLANK as usize, CORE_END.min(core / 2) as usize);
    let mut piece = [&seq[s - f..s + c], &seq[e - c..e + f]].concat();
    if minus {
        piece = dna::revcomp(&piece);
    }
    Some(OrientedCopy {
        locus_id: format!("{contig}:{start}-{end}"),
        left: f,
        right: f,
        seq: piece,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const OUT: &str = "   SW   perc perc perc  query      position in query           matching       repeat              position in repeat
score   div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID

  1000   10.0  0.0  0.0  chr1          101      600   (1000) + fam1           Unknown                 1    500     (0)      1
   400   20.0  0.0  0.0  chr1          801     1000    (600) C fam1           Unknown              (300)    200       1      2
   300   20.0  0.0  0.0  chr1         1101     1300    (300) C fam1           Unknown                (0)    500     301      2 *
   900    5.0  0.0  0.0  chr1         1401     5000      (0) + fam1           Unknown                 1    480     (20)      3
   200   15.0  0.0  0.0  chr2            1      100      (0) + fam2           LTR/Gypsy              10    109     (0)      4
";

    #[test]
    fn parses_both_strands_and_groups_insertions() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.out");
        fs::write(&p, OUT).unwrap();
        let fams = read_rm_out(&p).unwrap();
        assert_eq!(fams["fam1"].len(), 4);
        let c = &fams["fam1"][1];
        assert!(c.minus);
        assert_eq!(
            (c.start, c.end, c.rep_start, c.rep_end),
            (800, 1000, 1, 200)
        );
        let e = copy_evidence_from_rm(&fams["fam1"], 500, None, &CopyParams::default());
        // insertion 1: full; insertion 2: two fragments covering 1-200 and 301-500 (80%); 3: tandem-like span
        assert_eq!(e.loci, 3);
        assert_eq!(e.full_length, 2);
        assert!(e.mean_identity > 80.0 && e.mean_identity < 95.0);
    }

    #[test]
    fn cuts_flanks_and_calls_tsds_from_the_genome() {
        let mut x = 7u64;
        let mut rnd = |n: usize| {
            (0..n)
                .map(|_| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    b"ACGT"[(x >> 33) as usize % 4]
                })
                .collect::<Vec<u8>>()
        };
        let element = rnd(400);
        let mut chrom = Vec::new();
        let mut hits = Vec::new();
        for i in 0..6u64 {
            chrom.extend(rnd(300));
            let tsd = rnd(5);
            chrom.extend(&tsd);
            let start = chrom.len() as u64;
            if i % 2 == 0 {
                chrom.extend(&element);
            } else {
                chrom.extend(dna::revcomp(&element));
            }
            let end = chrom.len() as u64;
            chrom.extend(&tsd);
            hits.push(RmHit {
                contig: "chr1".into(),
                start,
                end,
                minus: i % 2 == 1,
                family: "f".into(),
                class: "Unknown".into(),
                rep_start: 1,
                rep_end: 400,
                divergence: 1.0,
                score: 1000,
                id: i,
            });
        }
        chrom.extend(rnd(300));
        let genome = HashMap::from([("chr1".to_owned(), chrom)]);
        let e = copy_evidence_from_rm(&hits, 400, Some(&genome), &CopyParams::default());
        assert_eq!((e.loci, e.full_length), (6, 6));
        let tsd = e.tsd.unwrap();
        assert_eq!(tsd.length, 5);
        assert_eq!(tsd.copies_evaluated, 6);
        assert!(e.flank_similarity.unwrap() < 0.5);
    }
}
