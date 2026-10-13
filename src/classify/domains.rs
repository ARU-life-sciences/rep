//! Protein-domain architecture of consensus ORFs (Pfam, HMMER `hmmscan`).
//!
//! A curated subset of Pfam models with known roles in transposable elements
//! is extracted once from a full `Pfam-A.hmm`. ORF hits are placed in coding
//! order along the dominant strand; the resulting architecture
//! (e.g. `GAG-PR-RT-RH-INT`) implies a class, and combinations that no known
//! class shows are reported as novel architectures.

use std::{
    collections::{BTreeMap, HashMap, HashSet},
    fs,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
    process::Command,
};

use serde::Serialize;

use crate::{Error, ErrorKind, Result};

/// Functional role of a TE-associated domain, written as a short code.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum Role {
    Gag,
    Protease,
    PaoProtease,
    ReverseTranscriptase,
    RnaseH,
    Integrase,
    Chromo,
    Endonuclease,
    GiyYig,
    TyrosineRecombinase,
    PolB,
    HelitronHelicase,
    /// Generic SF1 helicase (PIF1): in Helitrons but also host genes and other TEs.
    Helicase,
    Huh,
    /// DDE or other transposase of the named RepeatMasker superfamily.
    Transposase(&'static str),
    UlpProtease,
}

impl Role {
    pub fn code(self) -> String {
        match self {
            Self::Gag => "GAG".into(),
            Self::Protease => "PR".into(),
            Self::PaoProtease => "PR(Pao)".into(),
            Self::ReverseTranscriptase => "RT".into(),
            Self::RnaseH => "RH".into(),
            Self::Integrase => "INT".into(),
            Self::Chromo => "CHD".into(),
            Self::Endonuclease => "APE".into(),
            Self::GiyYig => "GIY".into(),
            Self::TyrosineRecombinase => "YR".into(),
            Self::PolB => "POLB".into(),
            Self::HelitronHelicase => "HEL".into(),
            Self::Helicase => "HELI".into(),
            Self::Huh => "HUH".into(),
            Self::Transposase(sf) => format!("TPase({sf})"),
            Self::UlpProtease => "ULP".into(),
        }
    }
}

/// Pfam models used, by name, with their TE role.
pub const TE_DOMAINS: &[(&str, Role)] = &[
    ("Retrotrans_gag", Role::Gag),
    ("Transposase_28", Role::Gag),
    ("Gag_p10", Role::Gag),
    ("Gag_p17", Role::Gag),
    ("Gag_p24", Role::Gag),
    ("Gag_p30", Role::Gag),
    ("Gag_MA", Role::Gag),
    ("zf-CCHC", Role::Gag),
    ("RVP", Role::Protease),
    ("RVP_2", Role::Protease),
    ("Asp_protease_2", Role::Protease),
    ("gag-asp_proteas", Role::Protease),
    ("Peptidase_A17", Role::PaoProtease),
    ("RVT_1", Role::ReverseTranscriptase),
    ("RVT_2", Role::ReverseTranscriptase),
    ("RVT_N", Role::ReverseTranscriptase),
    ("RVT_thumb", Role::ReverseTranscriptase),
    ("RNase_H", Role::RnaseH),
    ("RVT_3", Role::RnaseH),
    ("rve", Role::Integrase),
    ("Integrase_H2C2", Role::Integrase),
    ("zf-H2C2", Role::Integrase),
    ("gag_pre-integrs", Role::Integrase),
    ("Integrase_Zn", Role::Integrase),
    ("Chromo", Role::Chromo),
    ("Exo_endo_phos", Role::Endonuclease),
    ("Exo_endo_phos_2", Role::Endonuclease),
    ("GIY-YIG", Role::GiyYig),
    ("Phage_integrase", Role::TyrosineRecombinase),
    ("DNA_pol_B_2", Role::PolB),
    ("Helitron_like_N", Role::HelitronHelicase),
    ("PIF1", Role::Helicase),
    ("Rep_trans", Role::Huh),
    ("Peptidase_C48", Role::UlpProtease),
    ("DDE_1", Role::Transposase("TcMar-Pogo")),
    ("CENP-B_N", Role::Transposase("TcMar-Pogo")),
    ("HTH_Tnp_Tc5", Role::Transposase("TcMar-Pogo")),
    ("DDE_3", Role::Transposase("TcMar")),
    ("Transposase_1", Role::Transposase("TcMar")),
    ("HTH_48", Role::Transposase("TcMar")),
    ("HTH_SB", Role::Transposase("TcMar")),
    ("HTH_Tnp_Tc3_1", Role::Transposase("TcMar")),
    ("HTH_Tnp_Tc3_2", Role::Transposase("TcMar")),
    ("Dimer_Tnp_hAT", Role::Transposase("hAT")),
    ("hAT-like_RNase-H", Role::Transposase("hAT")),
    ("DBD_Tnp_Hermes", Role::Transposase("hAT")),
    ("Transposase_mut", Role::Transposase("MULE")),
    ("MULE", Role::Transposase("MULE")),
    ("DBD_Tnp_Mut", Role::Transposase("MULE")),
    ("Transposase_21", Role::Transposase("CMC-EnSpm")),
    ("Transposase_24", Role::Transposase("CMC-EnSpm")),
    ("DDE_Tnp_4", Role::Transposase("PIF-Harbinger")),
    ("Plant_tran", Role::Transposase("PIF-Harbinger")),
    ("DDE_Tnp_1_7", Role::Transposase("PiggyBac")),
    ("Tnp_P_element", Role::Transposase("P")),
    ("KDZ", Role::Transposase("Zisupton")),
    ("DDE_Tnp_IS1595", Role::Transposase("Merlin")),
    ("DDE_Tnp_1", Role::Transposase("DDE")),
];

pub fn role_of(domain: &str) -> Option<Role> {
    TE_DOMAINS
        .iter()
        .find(|(name, _)| *name == domain)
        .map(|&(_, role)| role)
}

/// Pfam families that occur in transposons as well as host genes (DNA-binding
/// motifs, zinc fingers, accessory domains of TE proteins). They neither mark
/// a TE nor a host gene, so the host-gene scan ignores them.
const AMBIGUOUS_PREFIXES: &[&str] = &[
    "zf-",
    "HTH_",
    "DDE_",
    "Transposase",
    "Tnp",
    "Dimer_Tnp",
    "DBD_Tnp",
    "RVT",
    "RVP",
    "Retrotrans",
    "rve",
    "Phage_",
    "DUF",
    "Herpes_",
    "UBN2",
    "Myb_",
    "CENP",
    "Chromo",
    "Exo_endo",
    "DNA_pol",
    "Helicase",
    "Peptidase_A",
    "RNase_H",
    "SWIM",
    "THAP",
    "MULE",
    "Gag",
    "gag",
    "Integrase",
    "Endonuclease",
    "PIF1",
    "AAA",
    "Rep_",
    "Ribonuc",
    "Mu-",
    "Tc5",
    "Tc3",
    "YqaJ",
];

/// Whether a Pfam family can belong to a TE (see [`AMBIGUOUS_PREFIXES`]).
pub fn is_te_associated(domain: &str) -> bool {
    role_of(domain).is_some() || AMBIGUOUS_PREFIXES.iter().any(|p| domain.starts_with(p))
}

/// Host-gene domains per family from an `hmmsearch --domtblout` of full Pfam
/// on consensus ORFs: Pfam families that are not TE-associated, with the
/// best domain score of each, strongest first.
pub fn parse_host_domains(path: &Path) -> Result<BTreeMap<String, Vec<(String, f64)>>> {
    let mut best: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();
    for line in BufReader::new(fs::File::open(path)?).lines() {
        let line = line?;
        let f = line.split_whitespace().collect::<Vec<_>>();
        if line.starts_with('#') || f.len() < 22 || is_te_associated(f[3]) {
            continue;
        }
        // hmmsearch: target is the ORF `<family>|orf<n>|<coords>`, query the domain
        let Some(family) = f[0].rsplitn(3, '|').nth(2) else {
            continue;
        };
        let score: f64 = f[13].parse()?;
        let e = best
            .entry(family.to_owned())
            .or_default()
            .entry(f[3].to_owned())
            .or_insert(score);
        *e = e.max(score);
    }
    Ok(best
        .into_iter()
        .map(|(family, domains)| {
            let mut v = domains.into_iter().collect::<Vec<_>>();
            v.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
            (family, v)
        })
        .collect())
}

/// Copy the TE domain models out of a full Pfam-A.hmm. Returns the names
/// missing from that Pfam release.
pub fn extract_te_models(pfam: &Path, out: &Path) -> Result<Vec<String>> {
    let wanted: HashSet<&str> = TE_DOMAINS.iter().map(|(n, _)| *n).collect();
    let mut found = HashSet::new();
    let mut writer = BufWriter::new(fs::File::create(out)?);
    let mut record = Vec::new();
    let mut keep = false;
    for line in BufReader::new(fs::File::open(pfam)?).lines() {
        let line = line?;
        if let Some(name) = line.strip_prefix("NAME  ") {
            keep = wanted.contains(name.trim());
            if keep {
                found.insert(name.trim().to_owned());
            }
        }
        record.push(line);
        if record.last().is_some_and(|l| l == "//") {
            if keep {
                for l in &record {
                    writeln!(writer, "{l}")?;
                }
            }
            record.clear();
            keep = false;
        }
    }
    writer.flush()?;
    let mut missing = wanted
        .into_iter()
        .filter(|n| !found.contains(*n))
        .map(str::to_owned)
        .collect::<Vec<_>>();
    missing.sort();
    if found.is_empty() {
        return Err(Error::new(ErrorKind::GenericCli(format!(
            "no TE domain models found in {}",
            pfam.display()
        ))));
    }
    Ok(missing)
}

pub fn hmmpress_command(exe: impl AsRef<Path>, hmm: &Path) -> Command {
    let mut c = Command::new(exe.as_ref());
    c.arg("-f").arg(hmm);
    c
}

/// `hmmscan` with Pfam gathering thresholds; per-domain hits to `domtblout`.
pub fn hmmscan_command(
    exe: impl AsRef<Path>,
    hmm: &Path,
    proteins: &Path,
    domtblout: &Path,
    threads: usize,
) -> Command {
    let mut c = Command::new(exe.as_ref());
    c.args(["--cut_ga", "--noali", "-o", "/dev/null", "--cpu"])
        .arg(threads.to_string())
        .arg("--domtblout")
        .arg(domtblout)
        .arg(hmm)
        .arg(proteins);
    c
}

/// One domain hit on a consensus ORF.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct DomainHit {
    pub family_id: String,
    pub domain: String,
    pub role: Role,
    pub orf: String,
    pub strand: char,
    pub score: f64,
    pub i_evalue: f64,
    /// Domain envelope in consensus nucleotides, 0-based half-open, measured
    /// along the coding strand (so order reflects translation order).
    pub coding_start: usize,
    pub coding_end: usize,
}

/// Parse an `hmmscan --domtblout` file of hits to ORFs named
/// `<family_id>|orf<n>|<strand>:<start>-<end>` (forward-strand nucleotide
/// coordinates); `lengths` gives each family's consensus length.
pub fn parse_domtblout(path: &Path, lengths: &HashMap<String, usize>) -> Result<Vec<DomainHit>> {
    let mut hits = Vec::new();
    for line in BufReader::new(fs::File::open(path)?).lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f = line.split_whitespace().collect::<Vec<_>>();
        if f.len() < 22 {
            return Err(Error::new(ErrorKind::GenericCli(format!(
                "malformed domtblout line: {line}"
            ))));
        }
        let Some(role) = role_of(f[0]) else { continue };
        let mut parts = f[3].rsplitn(3, '|');
        let (coords, orf_name, family) = (parts.next(), parts.next(), parts.next());
        let (Some(coords), Some(_), Some(family)) = (coords, orf_name, family) else {
            continue;
        };
        let (strand, range) = coords.split_at(1);
        let strand = strand.chars().next().unwrap_or('+');
        let (start, end) = range
            .trim_start_matches(':')
            .split_once('-')
            .map(|(s, e)| (s.parse::<usize>(), e.parse::<usize>()))
            .ok_or_else(|| Error::new(ErrorKind::GenericCli(format!("bad ORF name {}", f[3]))))?;
        let (start, end) = (start?, end?);
        let n = lengths.get(family).copied().unwrap_or(end);
        let env_from: usize = f[19].parse()?;
        let env_to: usize = f[20].parse()?;
        // ORF start along the coding strand
        let orf_coding_start = if strand == '+' { start } else { n - end };
        hits.push(DomainHit {
            family_id: family.to_owned(),
            domain: f[0].to_owned(),
            role,
            orf: f[3].to_owned(),
            strand,
            score: f[13].parse()?,
            i_evalue: f[12].parse()?,
            coding_start: orf_coding_start + 3 * (env_from - 1),
            coding_end: orf_coding_start + 3 * env_to,
        });
    }
    Ok(hits)
}

#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct DomainEvidence {
    /// Non-redundant hits on the dominant strand, in coding order.
    pub hits: Vec<DomainHit>,
    /// e.g. `GAG-PR-RT-RH-INT`
    pub architecture: String,
    /// Class implied by the architecture, if any.
    pub inferred: Option<String>,
    /// Why the architecture matches no known class, if it does not.
    pub novel: Option<String>,
}

/// Order, de-duplicate and interpret one family's domain hits.
pub fn interpret(hits: &[&DomainHit]) -> DomainEvidence {
    if hits.is_empty() {
        return DomainEvidence::default();
    }
    let mut strand_score: BTreeMap<char, f64> = BTreeMap::new();
    for h in hits {
        *strand_score.entry(h.strand).or_default() += h.score;
    }
    let strand = strand_score
        .iter()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(s, _)| *s)
        .expect("non-empty");
    // keep the best-scoring hit among those overlapping by more than half
    let mut candidates = hits
        .iter()
        .filter(|h| h.strand == strand)
        .copied()
        .collect::<Vec<_>>();
    candidates.sort_by(|a, b| b.score.total_cmp(&a.score));
    let mut kept: Vec<&DomainHit> = Vec::new();
    for h in candidates {
        let len = h.coding_end - h.coding_start;
        let redundant = kept.iter().any(|k| {
            let overlap = h
                .coding_end
                .min(k.coding_end)
                .saturating_sub(h.coding_start.max(k.coding_start));
            overlap * 2 > len.min(k.coding_end - k.coding_start)
        });
        if !redundant {
            kept.push(h);
        }
    }
    kept.sort_by_key(|h| h.coding_start);
    let mut codes: Vec<String> = Vec::new();
    for h in &kept {
        let code = h.role.code();
        if codes.last() != Some(&code) {
            codes.push(code);
        }
    }
    let (inferred, novel) = classify_architecture(&kept);
    DomainEvidence {
        hits: kept.into_iter().cloned().collect(),
        architecture: codes.join("-"),
        inferred,
        novel,
    }
}

/// Map an ordered set of domain hits to a class, or explain why it is novel.
fn classify_architecture(hits: &[&DomainHit]) -> (Option<String>, Option<String>) {
    let first = |role: Role| hits.iter().position(|h| h.role == role);
    let has = |role: Role| first(role).is_some();
    let transposases = hits
        .iter()
        .filter_map(|h| match h.role {
            Role::Transposase(sf) => Some(sf),
            _ => None,
        })
        .collect::<std::collections::BTreeSet<_>>();
    let rt = has(Role::ReverseTranscriptase);
    let yr = has(Role::TyrosineRecombinase);
    let polb = has(Role::PolB);
    let hel = has(Role::HelitronHelicase);

    // combinations that no described class shows
    let mut odd = Vec::new();
    if rt && !transposases.is_empty() {
        odd.push("RT+DDE-transposase");
    }
    if rt && polb {
        odd.push("RT+PolB");
    }
    if rt && hel {
        odd.push("RT+Helitron-helicase");
    }
    if yr && !transposases.is_empty() {
        odd.push("YR+DDE-transposase");
    }
    if yr && polb {
        odd.push("YR+PolB");
    }
    if transposases.len() > 1 && !transposases.contains("DDE") {
        odd.push("multiple-transposase-superfamilies");
    }
    if !odd.is_empty() {
        return (None, Some(odd.join(",")));
    }

    let inferred = if has(Role::PaoProtease) {
        // the Pao/BEL protease is diagnostic even when the RT is too divergent to detect
        Some("LTR/Pao")
    } else if rt {
        if yr {
            Some("LTR/DIRS")
        } else if has(Role::GiyYig) {
            Some("PLE/Penelope")
        } else if let (Some(int), Some(rt_pos)) =
            (first(Role::Integrase), first(Role::ReverseTranscriptase))
        {
            Some(if int < rt_pos {
                "LTR/Copia"
            } else {
                "LTR/Gypsy"
            })
        } else if has(Role::Endonuclease) {
            Some("LINE")
        } else {
            // RT without integrase/endonuclease: a retroelement of unclear order
            None
        }
    } else if polb {
        Some("DNA/Maverick")
    } else if hel || (has(Role::Huh) && has(Role::Helicase)) {
        Some("RC/Helitron")
    } else if yr {
        Some("DNA/Crypton")
    } else if !transposases.is_empty() {
        // prefer a superfamily-specific model over the generic DDE one
        return match transposases.iter().find(|sf| **sf != "DDE") {
            Some(sf) => (Some(format!("DNA/{sf}")), None),
            None => (Some("DNA".to_owned()), None),
        };
    } else {
        None
    };
    (inferred.map(str::to_owned), None)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(domain: &str, start: usize, score: f64) -> DomainHit {
        DomainHit {
            family_id: "f".into(),
            domain: domain.into(),
            role: role_of(domain).unwrap(),
            orf: "f|orf1|+:0-9000".into(),
            strand: '+',
            score,
            i_evalue: 1e-20,
            coding_start: start,
            coding_end: start + 300,
        }
    }

    fn arch(domains: &[(&str, usize)]) -> DomainEvidence {
        let hits = domains
            .iter()
            .map(|(d, s)| hit(d, *s, 50.0))
            .collect::<Vec<_>>();
        interpret(&hits.iter().collect::<Vec<_>>())
    }

    #[test]
    fn integrase_position_separates_gypsy_from_copia() {
        let gypsy = arch(&[
            ("Retrotrans_gag", 0),
            ("RVP", 600),
            ("RVT_1", 1200),
            ("RNase_H", 2000),
            ("rve", 2600),
        ]);
        assert_eq!(gypsy.architecture, "GAG-PR-RT-RH-INT");
        assert_eq!(gypsy.inferred.as_deref(), Some("LTR/Gypsy"));
        let copia = arch(&[
            ("Retrotrans_gag", 0),
            ("rve", 900),
            ("RVT_2", 1800),
            ("RNase_H", 2600),
        ]);
        assert_eq!(copia.inferred.as_deref(), Some("LTR/Copia"));
    }

    #[test]
    fn other_classes() {
        assert_eq!(
            arch(&[("RVT_1", 0), ("Peptidase_A17", 900), ("rve", 1800)])
                .inferred
                .as_deref(),
            Some("LTR/Pao")
        );
        assert_eq!(
            arch(&[("RVT_1", 0), ("Phage_integrase", 900)])
                .inferred
                .as_deref(),
            Some("LTR/DIRS")
        );
        assert_eq!(
            arch(&[("RVT_1", 0), ("GIY-YIG", 900)]).inferred.as_deref(),
            Some("PLE/Penelope")
        );
        assert_eq!(
            arch(&[("Exo_endo_phos", 0), ("RVT_1", 900)])
                .inferred
                .as_deref(),
            Some("LINE")
        );
        assert_eq!(
            arch(&[("DNA_pol_B_2", 0), ("rve", 900)])
                .inferred
                .as_deref(),
            Some("DNA/Maverick")
        );
        assert_eq!(
            arch(&[("Helitron_like_N", 0)]).inferred.as_deref(),
            Some("RC/Helitron")
        );
        assert_eq!(
            arch(&[("Phage_integrase", 0)]).inferred.as_deref(),
            Some("DNA/Crypton")
        );
        assert_eq!(
            arch(&[("Dimer_Tnp_hAT", 0)]).inferred.as_deref(),
            Some("DNA/hAT")
        );
        assert_eq!(arch(&[("DDE_Tnp_1", 0)]).inferred.as_deref(), Some("DNA"));
        assert_eq!(arch(&[("RVT_1", 0)]).inferred, None);
        assert_eq!(
            arch(&[("Peptidase_A17", 0), ("rve", 900)])
                .inferred
                .as_deref(),
            Some("LTR/Pao")
        );
        // a generic PIF1 helicase is not a Helitron on its own
        assert_eq!(arch(&[("PIF1", 0)]).inferred, None);
        assert_eq!(
            arch(&[("Tnp_P_element", 0), ("PIF1", 900)])
                .inferred
                .as_deref(),
            Some("DNA/P")
        );
        assert_eq!(
            arch(&[("Rep_trans", 0), ("PIF1", 900)]).inferred.as_deref(),
            Some("RC/Helitron")
        );
    }

    #[test]
    fn unprecedented_combinations_are_novel() {
        let e = arch(&[("RVT_1", 0), ("Dimer_Tnp_hAT", 900)]);
        assert_eq!(
            (e.inferred, e.novel.as_deref()),
            (None, Some("RT+DDE-transposase"))
        );
        let e = arch(&[("Phage_integrase", 0), ("DNA_pol_B_2", 900)]);
        assert_eq!(e.novel.as_deref(), Some("YR+PolB"));
    }

    #[test]
    fn overlapping_hits_keep_the_best_and_minus_strand_is_ignored_when_minor() {
        let mut weak = hit("RVT_2", 10, 20.0);
        weak.coding_end = 300;
        let mut other_strand = hit("Dimer_Tnp_hAT", 5000, 10.0);
        other_strand.strand = '-';
        let hits = [
            hit("RVT_1", 0, 80.0),
            weak,
            other_strand,
            hit("rve", 900, 60.0),
        ];
        let e = interpret(&hits.iter().collect::<Vec<_>>());
        assert_eq!(e.architecture, "RT-INT");
        assert_eq!(e.hits[0].domain, "RVT_1");
        assert!(e.novel.is_none());
    }

    #[test]
    fn parses_domtblout_and_maps_minus_strand_orfs_to_coding_order() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("d.tbl");
        // ORF on the minus strand spanning consensus 100..1000 of a 2000 bp family
        let row = |domain: &str, from: usize, to: usize| {
            format!("{domain} PF0 200 fam#Unknown|orf1|-:100-1000 - 300 1e-30 100.0 0.1 1 1 1e-30 1e-30 99.0 0.1 1 200 {from} {to} {from} {to} 0.9 desc words\n")
        };
        fs::write(
            &path,
            format!("# header\n{}{}", row("RVT_1", 1, 100), row("rve", 150, 250)),
        )
        .unwrap();
        let lengths = HashMap::from([("fam#Unknown".to_owned(), 2000usize)]);
        let hits = parse_domtblout(&path, &lengths).unwrap();
        assert_eq!(hits.len(), 2);
        // coding strand starts at 2000 - 1000 = 1000 in reverse-complement coordinates
        assert_eq!((hits[0].coding_start, hits[0].coding_end), (1000, 1300));
        assert_eq!(hits[1].coding_start, 1000 + 3 * 149);
        assert_eq!(hits[0].family_id, "fam#Unknown");
    }

    #[test]
    fn host_domains_skip_te_associated_families() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("host.domtblout");
        let row = |orf: &str, dom: &str, score: f64| {
            format!("{orf} - 300 {dom} PF0 200 1e-30 {score} 0.1 1 1 1e-30 1e-30 {score} 0.1 1 200 1 200 1 200 0.9 -\n")
        };
        fs::write(
            &p,
            [
                "# comment\n".to_owned(),
                row("f#Unknown|orf1|+:0-900", "Actin", 250.0),
                row("f#Unknown|orf2|-:0-900", "Actin", 300.0),
                row("f#Unknown|orf1|+:0-900", "zf-C2H2", 40.0),
                row("f#Unknown|orf1|+:0-900", "RVT_1", 90.0),
                row("g#Unknown|orf1|+:0-900", "7tm_1", 80.0),
            ]
            .concat(),
        )
        .unwrap();
        let host = parse_host_domains(&p).unwrap();
        assert_eq!(host["f#Unknown"], vec![("Actin".to_owned(), 300.0)]);
        assert_eq!(host["g#Unknown"], vec![("7tm_1".to_owned(), 80.0)]);
        assert!(is_te_associated("DDE_Tnp_1") && !is_te_associated("Ribosomal_S2"));
    }
}
