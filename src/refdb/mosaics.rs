//! Reference consensi that carry a copy of an element of another class.
//!
//! Some reference consensi contain an embedded element of a different class,
//! e.g. a Helitron inside a Repbase Kolobok consensus, presumably inserted in
//! the copies the consensus was built from. A family homologous to the
//! embedded element then matches the host reference and inherits its label.
//! The nucleotide library is searched against itself. A stretch of reference
//! R (at most half of it) is a foreign segment when more references of
//! another class than of R's own class match it locally; `rep classify` then
//! attributes hits within it to the embedded class. The vote keeps pieces
//! that classes legitimately share (e.g. the LINE-derived tails of SINEs)
//! with their host, as relatives of the host match there too.

use std::{
    collections::{BTreeMap, HashMap},
    fs,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    process::Command,
};

use bio::io::fasta;
use serde::{Deserialize, Serialize};

use super::{cli_error, BuildInfo, RefDb};
use crate::{
    classify::{canonical_class, is_non_te, split_label},
    command_runner::run_checked,
    copy_recovery::executable,
    CommandRunner, Result,
};

#[derive(Clone, Debug, Serialize)]
pub struct MosaicOptions {
    /// Fraction of the embedded reference the alignment must cover (0: any
    /// local match counts; embedded copies are usually fragments).
    pub min_contained: f64,
    /// Largest fraction of the host reference a foreign segment may cover
    /// (beyond it, the two references are simply similar).
    pub max_host_fraction: f64,
    pub min_length: u64,
    pub threads: usize,
}

impl Default for MosaicOptions {
    fn default() -> Self {
        Self {
            min_contained: 0.0,
            max_host_fraction: 0.5,
            min_length: 80,
            threads: 1,
        }
    }
}

/// One foreign segment of a reference (1-based inclusive coordinates).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Mosaic {
    pub reference: String,
    pub start: u64,
    pub end: u64,
    /// Label of the embedded element (most common among the contained references).
    pub class: String,
    /// Distinct references of that class contained in the segment.
    pub support: usize,
    pub example: String,
    /// Contained references (not written out).
    #[serde(skip)]
    pub guests: Vec<String>,
}

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct MosaicInfo {
    pub references: usize,
    pub mosaic_references: usize,
    pub segments: usize,
    pub by_pair: BTreeMap<String, usize>,
    pub source_fingerprint: String,
}

impl RefDb {
    pub fn mosaics(&self) -> PathBuf {
        self.build_dir().join("mosaics.tsv")
    }
}

/// One alignment of the self-search: reference `host` contains `guest`.
struct Containment {
    guest: String,
    guest_class: String,
    start: u64,
    end: u64,
    host_length: u64,
}

/// SINE 3' tails derive from LINE 3' ends, so the two classes share that
/// stretch legitimately: a LINE segment in a SINE's 3' part, or a SINE
/// segment at a LINE's 3' end, is not an embedded element.
fn shared_sine_line_end(host_class: &str, guest_class: &str, end: u64, host_length: u64) -> bool {
    let tail = |fraction: f64| end as f64 >= fraction * host_length as f64;
    match (host_class, guest_class) {
        ("SINE", "LINE") => tail(0.6),
        ("LINE", "SINE") => tail(0.9),
        _ => false,
    }
}

/// Search the nucleotide library against itself and write `build/mosaics.tsv`.
pub fn find_mosaics(
    db: &RefDb,
    options: &MosaicOptions,
    runner: &dyn CommandRunner,
) -> Result<MosaicInfo> {
    let library = db.nucleotide_library();
    let work = db.build_dir().join("mosaics_work");
    if work.exists() {
        fs::remove_dir_all(&work)?;
    }
    fs::create_dir_all(&work)?;
    let hits = work.join("self.m8");
    let mut search = Command::new(executable("mmseqs"));
    search
        .arg("easy-search")
        .arg(&library)
        .arg(&library)
        .arg(&hits)
        .arg(work.join("tmp"))
        .args([
            "--search-type",
            "3",
            "-e",
            "1e-10",
            "--format-output",
            "query,target,qstart,qend,qlen,tstart,tend,tlen,pident,bits",
            "-v",
            "1",
            "--threads",
            &options.threads.to_string(),
        ]);
    run_checked(runner, &mut search, "mmseqs easy-search")?;
    let mut info = MosaicInfo {
        references: fasta::Reader::from_file(&library)?.records().count(),
        ..Default::default()
    };
    let mosaics = mosaics_from_hits(&hits, options)?;
    write_mosaics(&db.mosaics(), &mosaics)?;
    info.mosaic_references = mosaics
        .iter()
        .map(|m| m.reference.as_str())
        .collect::<std::collections::BTreeSet<_>>()
        .len();
    info.segments = mosaics.len();
    for m in &mosaics {
        let host = class_of(&m.reference).unwrap_or_default();
        *info
            .by_pair
            .entry(format!("{} in {host}", canonical_class(&m.class)))
            .or_default() += 1;
    }
    info.source_fingerprint = fs::read(db.build_dir().join("build.json"))
        .ok()
        .and_then(|b| serde_json::from_slice::<BuildInfo>(&b).ok())
        .map(|b| b.fingerprint)
        .unwrap_or_default();
    fs::write(
        db.build_dir().join("mosaics.json"),
        serde_json::to_vec_pretty(&info).map_err(|e| cli_error(e.to_string()))?,
    )?;
    fs::remove_dir_all(&work).ok();
    Ok(info)
}

/// Top-level TE class of a `name#label` reference, if it has one.
fn class_of(id: &str) -> Option<String> {
    let (_, label) = split_label(id);
    let label = label?;
    let class = canonical_class(&label);
    (!class.is_empty() && class != "Unknown" && !is_non_te(&label)).then(|| class.to_owned())
}

/// Foreign segments from a self-search table (`query,target,qstart,qend,qlen,
/// tstart,tend,tlen,pident,bits`), where the query is the host.
pub fn mosaics_from_hits(path: &Path, o: &MosaicOptions) -> Result<Vec<Mosaic>> {
    let mut contained: BTreeMap<String, Vec<Containment>> = BTreeMap::new();
    for line in BufReader::new(fs::File::open(path)?).lines() {
        let line = line?;
        let f = line.split('\t').collect::<Vec<_>>();
        if f.len() < 10 || f[0] == f[1] {
            continue;
        }
        let num = |i: usize| f[i].parse::<u64>().ok();
        let (Some(qs), Some(qe), Some(qlen), Some(ts), Some(te), Some(tlen)) =
            (num(2), num(3), num(4), num(5), num(6), num(7))
        else {
            continue;
        };
        let (start, end) = (qs.min(qe), qs.max(qe));
        let span = end - start + 1;
        let guest_span = ts.abs_diff(te) + 1;
        if span < o.min_length
            || (guest_span as f64) < o.min_contained * tlen as f64
            || span as f64 > o.max_host_fraction * qlen as f64
        {
            continue;
        }
        let Some(guest_class) = class_of(f[1]) else {
            continue;
        };
        contained
            .entry(f[0].to_owned())
            .or_default()
            .push(Containment {
                guest: f[1].to_owned(),
                guest_class,
                start,
                end,
                host_length: qlen,
            });
    }

    let mut mosaics = Vec::new();
    for (host, all) in &contained {
        let Some(host_class) = class_of(host) else {
            continue;
        };
        // merge overlapping containments per class into segments
        let mut by_class: HashMap<&str, Vec<&Containment>> = HashMap::new();
        for g in all {
            by_class.entry(&g.guest_class).or_default().push(g);
        }
        for (class, mut list) in by_class {
            if class == host_class {
                continue;
            }
            list.sort_by_key(|g| (g.start, g.end));
            let mut segments: Vec<(u64, u64, Vec<&Containment>)> = Vec::new();
            for g in list {
                match segments.last_mut() {
                    Some(s) if g.start <= s.1 => {
                        s.1 = s.1.max(g.end);
                        s.2.push(g);
                    }
                    _ => segments.push((g.start, g.end, vec![g])),
                }
            }
            for (start, end, members) in segments {
                if shared_sine_line_end(&host_class, class, end, members[0].host_length) {
                    continue;
                }
                let mut guests = members.iter().map(|g| g.guest.as_str()).collect::<Vec<_>>();
                guests.sort_unstable();
                guests.dedup();
                // references of the host's own class contained in the same
                // stretch argue that it belongs to the host
                let mut own = all
                    .iter()
                    .filter(|g| {
                        // overlapping at least half of the segment or of the relative's match
                        let overlap = (g.end.min(end) + 1).saturating_sub(g.start.max(start));
                        g.guest_class == host_class
                            && 2 * overlap >= (end - start + 1).min(g.end - g.start + 1)
                    })
                    .map(|g| g.guest.as_str())
                    .collect::<Vec<_>>();
                own.sort_unstable();
                own.dedup();
                if guests.len() > own.len() {
                    // the most common full label among the embedded references
                    let mut labels: BTreeMap<String, usize> = BTreeMap::new();
                    for g in &guests {
                        if let (_, Some(l)) = split_label(g) {
                            *labels.entry(l).or_default() += 1;
                        }
                    }
                    let label = labels
                        .into_iter()
                        .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(&a.0)))
                        .map_or_else(|| class.to_owned(), |(l, _)| l);
                    mosaics.push(Mosaic {
                        reference: host.clone(),
                        start,
                        end,
                        class: label,
                        support: guests.len(),
                        example: guests[0].to_owned(),
                        guests: guests.iter().map(|g| (*g).to_owned()).collect(),
                    });
                }
            }
        }
    }
    // two references sharing a stretch can each look like the other's host:
    // keep the better-supported direction, neither on a tie
    let support_of = |host: &str, guest: &str| {
        mosaics
            .iter()
            .filter(|m| m.reference == host && m.guests.iter().any(|g| g == guest))
            .map(|m| m.support)
            .max()
    };
    let keep = mosaics
        .iter()
        .map(|m| {
            m.guests
                .iter()
                .all(|g| support_of(g, &m.reference).is_none_or(|reverse| m.support > reverse))
        })
        .collect::<Vec<_>>();
    let mut mosaics = mosaics
        .into_iter()
        .zip(keep)
        .filter_map(|(m, k)| k.then_some(m))
        .collect::<Vec<_>>();
    mosaics.sort_by(|a, b| (&a.reference, a.start).cmp(&(&b.reference, b.start)));
    Ok(mosaics)
}

fn write_mosaics(path: &Path, mosaics: &[Mosaic]) -> Result<()> {
    let mut text = String::from("reference\tstart\tend\tclass\tsupport\texample\n");
    for m in mosaics {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            m.reference, m.start, m.end, m.class, m.support, m.example
        ));
    }
    let tmp = path.with_extension("tsv.tmp");
    fs::write(&tmp, text)?;
    fs::rename(tmp, path)?;
    Ok(())
}

/// Foreign segments by reference, as written by [`find_mosaics`].
pub fn read_mosaics(path: &Path) -> Result<HashMap<String, Vec<Mosaic>>> {
    let mut by_reference: HashMap<String, Vec<Mosaic>> = HashMap::new();
    for line in BufReader::new(fs::File::open(path)?).lines().skip(1) {
        let line = line?;
        let f = line.split('\t').collect::<Vec<_>>();
        if f.len() < 6 {
            continue;
        }
        let (Ok(start), Ok(end), Ok(support)) = (f[1].parse(), f[2].parse(), f[4].parse()) else {
            continue;
        };
        by_reference
            .entry(f[0].to_owned())
            .or_default()
            .push(Mosaic {
                reference: f[0].to_owned(),
                start,
                end,
                class: f[3].to_owned(),
                support,
                example: f[5].to_owned(),
                guests: Vec::new(),
            });
    }
    Ok(by_reference)
}

/// The embedded class of a reference stretch (1-based inclusive, either
/// order), when at least half of it lies in one foreign segment.
pub fn foreign_class<'a>(
    mosaics: &'a HashMap<String, Vec<Mosaic>>,
    reference: &str,
    start: u64,
    end: u64,
) -> Option<&'a Mosaic> {
    let (s, e) = (start.min(end), start.max(end));
    mosaics.get(reference)?.iter().find(|m| {
        let overlap = (e.min(m.end) + 1).saturating_sub(s.max(m.start));
        2 * overlap >= e - s + 1
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(q: &str, t: &str, qs: u64, qe: u64, qlen: u64, ts: u64, te: u64, tlen: u64) -> String {
        format!("{q}\t{t}\t{qs}\t{qe}\t{qlen}\t{ts}\t{te}\t{tlen}\t85.0\t200\n")
    }

    #[test]
    fn embedded_elements_of_other_classes_are_foreign_segments() {
        let dir = tempfile::tempdir().unwrap();
        let hits = dir.path().join("self.m8");
        let k = "KolobokE-1#DNA/Kolobok";
        let rows = [
            // two Helitrons lie inside the Kolobok consensus at 5363-5559
            row(k, "Hel-1#RC/Helitron", 5363, 5559, 5826, 1, 197, 263),
            row(k, "Hel-2#RC/Helitron", 5400, 5600, 5826, 10, 200, 250),
            // a fragment of a long Helitron, outvoted by two Kolobok relatives
            row(k, "Hel-3#RC/Helitron", 100, 300, 5826, 1, 200, 2000),
            row(k, "KolobokN-1#DNA/Kolobok", 90, 310, 5826, 1, 220, 230),
            row(k, "KolobokN-2#DNA/Kolobok", 100, 300, 5826, 1, 200, 210),
            // the reverse direction: the Helitron is not a mosaic of Kolobok
            row("Hel-1#RC/Helitron", k, 1, 197, 263, 5363, 5559, 5826),
            // a SINE tail at a LINE's 3' end, where other LINEs match too
            row(
                "L2-1#LINE/L2",
                "MIR#SINE/MIR",
                3000,
                3100,
                3300,
                160,
                260,
                260,
            ),
            row(
                "L2-1#LINE/L2",
                "L2-2#LINE/L2",
                2800,
                3300,
                3300,
                1,
                500,
                3200,
            ),
            // same class, and Unknown guests, are ignored
            row(k, "KolobokN#DNA/Kolobok", 2000, 2200, 5826, 1, 200, 210),
            row(k, "x#Unknown", 3000, 3200, 5826, 1, 200, 200),
        ]
        .concat();
        fs::write(&hits, rows).unwrap();
        let m = mosaics_from_hits(&hits, &MosaicOptions::default()).unwrap();
        assert_eq!(m.len(), 1, "{m:?}");
        assert_eq!(
            (
                m[0].reference.as_str(),
                m[0].start,
                m[0].end,
                m[0].class.as_str(),
                m[0].support
            ),
            (k, 5363, 5600, "RC/Helitron", 2)
        );
        let path = dir.path().join("mosaics.tsv");
        write_mosaics(&path, &m).unwrap();
        let read = read_mosaics(&path).unwrap();
        let strip = |v: &[Mosaic]| {
            v.iter()
                .map(|m| Mosaic {
                    guests: Vec::new(),
                    ..m.clone()
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(read[k], strip(&m));
        assert_eq!(
            foreign_class(&read, k, 5559, 5380).unwrap().class,
            "RC/Helitron"
        );
        assert!(foreign_class(&read, k, 1, 5826).is_none());
        assert!(foreign_class(&read, "other", 1, 10).is_none());
    }

    #[test]
    fn own_class_guests_outvote_a_foreign_one() {
        let dir = tempfile::tempdir().unwrap();
        let hits = dir.path().join("self.m8");
        let h = "Mariner-1#DNA/TcMar";
        let rows = [
            row(h, "Hel-1#RC/Helitron", 100, 300, 1300, 1, 200, 210),
            row(h, "nMar-1#DNA/TcMar", 100, 300, 1300, 1, 200, 200),
            row(h, "nMar-2#DNA/TcMar", 110, 290, 1300, 1, 180, 190),
        ]
        .concat();
        fs::write(&hits, rows).unwrap();
        assert!(mosaics_from_hits(&hits, &MosaicOptions::default())
            .unwrap()
            .is_empty());
    }

    #[test]
    fn mutual_mosaics_keep_the_better_supported_direction() {
        let dir = tempfile::tempdir().unwrap();
        let hits = dir.path().join("self.m8");
        let (p, h1, h2) = ("P-2#DNA/P", "Hel-1#RC/Helitron", "Hel-2#RC/Helitron");
        let rows = [
            // the P consensus carries a Helitron fragment matched by two Helitrons;
            // each Helitron's matching stretch is matched by one P
            row(p, h1, 5228, 5388, 14408, 300, 460, 1150),
            row(p, h2, 5230, 5380, 14408, 400, 550, 700),
            row(h1, p, 300, 460, 1150, 5228, 5388, 14408),
            row(h2, p, 400, 550, 700, 5230, 5380, 14408),
        ]
        .concat();
        fs::write(&hits, rows).unwrap();
        let m = mosaics_from_hits(&hits, &MosaicOptions::default()).unwrap();
        assert_eq!(m.len(), 1, "{m:?}");
        assert_eq!(
            (m[0].reference.as_str(), m[0].class.as_str()),
            (p, "RC/Helitron")
        );
        // one against one: ambiguous, so neither
        fs::write(
            &hits,
            [
                row(p, h1, 5228, 5388, 14408, 300, 460, 1150),
                row(h1, p, 300, 460, 1150, 5228, 5388, 14408),
            ]
            .concat(),
        )
        .unwrap();
        assert!(mosaics_from_hits(&hits, &MosaicOptions::default())
            .unwrap()
            .is_empty());
    }

    #[test]
    fn shared_sine_line_ends_are_not_mosaics() {
        let dir = tempfile::tempdir().unwrap();
        let hits = dir.path().join("self.m8");
        let rows = [
            // the LINE-derived 3' tail of a SINE, and a LINE 3' end matched by SINE tails
            row(
                "SINE2-1#SINE/tRNA",
                "BovB#LINE/RTE",
                124,
                232,
                240,
                3000,
                3110,
                3200,
            ),
            row(
                "L2-1#LINE/L2",
                "MIR#SINE/MIR",
                3150,
                3290,
                3300,
                120,
                260,
                260,
            ),
            // a SINE inserted in the middle of a LINE consensus is a mosaic
            row("L1-1#LINE/L1", "B2#SINE/B2", 1288, 1484, 6000, 1, 190, 200),
        ]
        .concat();
        fs::write(&hits, rows).unwrap();
        let m = mosaics_from_hits(&hits, &MosaicOptions::default()).unwrap();
        assert_eq!(m.len(), 1, "{m:?}");
        assert_eq!(m[0].reference, "L1-1#LINE/L1");
    }
}
