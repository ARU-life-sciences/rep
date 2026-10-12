//! Combine homology and structural evidence into a per-family disposition.
//!
//! The rules are deliberately transparent: every disposition lists the
//! evidence that produced it, and anything not explained is `unresolved` with
//! the reason recorded.

use serde::Serialize;

use super::{
    copies::{expected_tsd, superfamilies_with_tsd, CopyEvidence},
    domains::DomainEvidence,
    homology::{ClassSupport, HomologySummary},
    structure::StructuralEvidence,
    taxonomy::top_level,
};

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum DispositionKind {
    KnownTe,
    Chimeric,
    SimpleRepeat,
    Satellite,
    LowComplexity,
    /// Matches a non-TE reference: RNA genes, multicopy host genes, artefacts.
    HostOrArtefact,
    /// A confirmed dispersed, multicopy repeat with no TE hallmark yet.
    InterspersedRepeat,
    NovelCandidate,
    Unresolved,
}

impl DispositionKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::KnownTe => "known_te",
            Self::Chimeric => "chimeric",
            Self::SimpleRepeat => "simple_repeat",
            Self::Satellite => "satellite",
            Self::LowComplexity => "low_complexity",
            Self::HostOrArtefact => "host_or_artefact",
            Self::InterspersedRepeat => "interspersed_repeat",
            Self::NovelCandidate => "novel_candidate",
            Self::Unresolved => "unresolved",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Confidence {
    Low,
    Medium,
    High,
}

impl Confidence {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Low => "low",
            Self::Medium => "medium",
            Self::High => "high",
        }
    }

    fn lower(self) -> Self {
        match self {
            Self::High => Self::Medium,
            _ => Self::Low,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct DecisionThresholds {
    /// Minimum bitscore for a protein hit to count as homology.
    pub min_bitscore: f64,
    /// Bitscore and consensus coverage needed for high-confidence homology.
    pub strong_bitscore: f64,
    pub strong_coverage: f64,
    /// Tandem periods at or below this are simple repeats, above are satellites.
    pub simple_repeat_max_period: usize,
    /// Minimum number of tandem units for a tandem call.
    pub min_tandem_copies: f64,
    /// ORF length (aa) regarded as substantial coding capacity.
    pub coding_orf_aa: usize,
    /// Non-coding poly-A families shorter than this are SINE-like.
    pub sine_max_length: usize,
    /// Minimum bitscore and consensus coverage for nucleotide homology.
    pub nucl_min_bitscore: f64,
    pub nucl_min_coverage: f64,
    /// Coverage for high-confidence nucleotide homology.
    pub nucl_strong_coverage: f64,
}

impl Default for DecisionThresholds {
    fn default() -> Self {
        Self {
            min_bitscore: 50.0,
            strong_bitscore: 200.0,
            strong_coverage: 0.2,
            simple_repeat_max_period: 10,
            min_tandem_copies: 3.0,
            coding_orf_aa: 300,
            sine_max_length: 700,
            nucl_min_bitscore: 50.0,
            nucl_min_coverage: 0.1,
            nucl_strong_coverage: 0.6,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct Disposition {
    pub kind: DispositionKind,
    /// RepeatMasker label to write into the reclassified library.
    pub label: String,
    pub confidence: Confidence,
    pub evidence: Vec<String>,
    pub warnings: Vec<String>,
}

/// Structural hallmarks expected for each top-level class. A class with no
/// entry imposes no expectation.
fn structure_consistent(class: &str, s: &StructuralEvidence) -> Option<bool> {
    match top_level(class) {
        "LTR" if class.contains("DIRS") || class.contains("Ngaro") => None,
        "LTR" => Some(s.ltr.is_some() || s.tg_ca_termini),
        "DNA" if class.contains("Maverick") || class.contains("Crypton") => None,
        "DNA" => Some(s.tir.is_some()),
        "LINE" => Some(s.tail.is_some() || s.ltr.is_none() && s.tir.is_none()),
        _ => None,
    }
}

/// Labels that are repeats or reference sequences but not transposable elements.
pub fn is_non_te(label: &str) -> bool {
    matches!(
        top_level(label),
        "Satellite"
            | "Simple_repeat"
            | "Low_complexity"
            | "rRNA"
            | "tRNA"
            | "snRNA"
            | "scRNA"
            | "srpRNA"
            | "RNA"
            | "Other"
            | "ARTEFACT"
    )
}

/// A reference family that is itself unclassified beyond its order, such as
/// our own earlier `LTR/Unknown` novel candidates.
fn is_open_label(label: &str) -> bool {
    label.ends_with("/Unknown")
}

fn describe(source: &str, best: &ClassSupport) -> String {
    format!(
        "{source}:{}:{}(bits={:.0},cov={:.2},id={:.1})",
        best.class, best.best_subject, best.best_bitscore, best.query_coverage, best.best_identity
    )
}

pub fn decide(
    homology: &HomologySummary,
    nucleotide: &HomologySummary,
    domains: &DomainEvidence,
    structure: &StructuralEvidence,
    t: &DecisionThresholds,
) -> Disposition {
    let mut evidence = structural_notes(structure, t);
    let mut warnings = Vec::new();
    if !domains.architecture.is_empty() {
        evidence.push(format!("domains:{}", domains.architecture));
    }
    if let Some(reason) = &domains.novel {
        warnings.push(format!("novel_domain_architecture:{reason}"));
    }
    let protein_best = homology.classes.iter().find(|c| !is_non_te(&c.class));
    let nucl_best = nucleotide
        .best()
        .filter(|b| b.query_coverage >= t.nucl_min_coverage);
    if let Some(b) = nucl_best {
        evidence.insert(0, describe("nucleotide", b));
    }

    // 0. Nucleotide match to a non-TE reference, unless TE proteins say otherwise
    if let Some(b) = nucl_best.filter(|b| is_non_te(&b.class) && protein_best.is_none()) {
        let kind = match top_level(&b.class) {
            "Satellite" => DispositionKind::Satellite,
            "Simple_repeat" => DispositionKind::SimpleRepeat,
            "Low_complexity" => DispositionKind::LowComplexity,
            _ => DispositionKind::HostOrArtefact,
        };
        return Disposition {
            kind,
            label: b.class.clone(),
            confidence: if b.query_coverage >= t.nucl_strong_coverage {
                Confidence::High
            } else {
                Confidence::Medium
            },
            evidence,
            warnings,
        };
    }

    // 1. Homology to a classified TE: proteins first, then nucleotide consensi
    let te_nucl = nucl_best.filter(|b| !is_non_te(&b.class));
    let homologous = match (protein_best, te_nucl) {
        (Some(p), _) => Some((p, "protein", homology, t.strong_coverage)),
        (None, Some(n)) => Some((n, "nucleotide", nucleotide, t.nucl_strong_coverage)),
        (None, None) => None,
    };
    if let Some((best, source, summary, strong_coverage)) = homologous {
        if source == "protein" {
            evidence.insert(0, describe("protein", best));
        }
        let mut confidence =
            if best.best_bitscore >= t.strong_bitscore && best.query_coverage >= strong_coverage {
                Confidence::High
            } else {
                Confidence::Medium
            };
        if let (Some(p), Some(n)) = (protein_best, te_nucl) {
            if top_level(&p.class) == top_level(&n.class) {
                evidence.push(format!("nucleotide_agrees:{}", n.class));
            } else if n.query_coverage >= t.nucl_strong_coverage {
                warnings.push(format!("nucleotide_class_conflict:{}", n.class));
                confidence = confidence.lower();
            }
        }
        if source == "nucleotide" && (best.best_bitscore < 100.0 || best.query_coverage < 0.3) {
            // weak nucleotide evidence: accepted, but only as a low-confidence call
            warnings.push("weak_nucleotide_homology".into());
            confidence = Confidence::Low;
        }
        if source == "nucleotide" && structure.length < 500 && best.best_subject_coverage < 0.1 {
            // a short family matching a small piece of a long reference
            warnings.push(format!(
                "short_match_within_long_reference:{:.2}",
                best.best_subject_coverage
            ));
            confidence = Confidence::Low;
        }
        if let Some(inferred) = &domains.inferred {
            if top_level(inferred) == top_level(&best.class) {
                evidence.push(format!("domains_agree:{inferred}"));
            } else {
                warnings.push(format!("domain_class_conflict:{inferred}"));
                confidence = confidence.lower();
            }
        }
        match structure_consistent(&best.class, structure) {
            Some(true) => evidence.push(format!("structure_consistent_with:{}", best.class)),
            Some(false) => {
                warnings.push(format!("expected_structure_missing_for:{}", best.class));
                confidence = confidence.lower();
            }
            None => {}
        }
        if structure.tandem.is_some() {
            warnings.push("tandem_periodicity_in_te_consensus".into());
        }
        if !summary.chimeric_classes.is_empty() {
            warnings.push(format!(
                "{source}_homology_to_multiple_classes:{}",
                summary.chimeric_classes.join(",")
            ));
            // SINEs (tRNA head + LINE tail) and Helitrons (captured fragments)
            // are composite by nature: hits to several classes are expected, and
            // the composite's own class is the right call
            for composite in ["SINE", "RC"] {
                if summary.chimeric_classes.iter().any(|c| c == composite) {
                    if let Some(own) = summary
                        .classes
                        .iter()
                        .find(|c| top_level(&c.class) == composite)
                    {
                        warnings.push(format!("composite_element:{}", own.class));
                        return Disposition {
                            kind: DispositionKind::KnownTe,
                            label: own.class.clone(),
                            confidence: Confidence::Low,
                            evidence,
                            warnings,
                        };
                    }
                }
            }
            return Disposition {
                kind: DispositionKind::Chimeric,
                label: best.class.clone(),
                confidence: Confidence::Low,
                evidence,
                warnings,
            };
        }
        if is_open_label(&best.class) {
            // a recurrence of a family that is itself only known to order level
            evidence.push(format!(
                "recurs_with_reference_family:{}",
                best.best_subject
            ));
            return Disposition {
                kind: DispositionKind::NovelCandidate,
                label: best.class.clone(),
                confidence: confidence.lower(),
                evidence,
                warnings,
            };
        }
        return Disposition {
            kind: DispositionKind::KnownTe,
            label: best.class.clone(),
            confidence,
            evidence,
            warnings,
        };
    }
    if let Some(subject) = &homology.best_unclassified_subject {
        evidence.insert(
            0,
            format!(
                "protein_unclassified_subject:{subject}(bits={:.0})",
                homology.best_unclassified_bitscore
            ),
        );
        warnings.push("homology_only_to_unclassified_te_proteins".into());
    }

    // 2. Protein-domain architecture without homology to known families
    if let Some(reason) = &domains.novel {
        return Disposition {
            kind: DispositionKind::NovelCandidate,
            label: "Unknown".into(),
            confidence: Confidence::Medium,
            evidence: {
                evidence.insert(0, format!("novel_domain_architecture:{reason}"));
                evidence
            },
            warnings,
        };
    }
    if let Some(inferred) = &domains.inferred {
        let mut confidence = Confidence::Medium;
        match structure_consistent(inferred, structure) {
            Some(true) => evidence.push(format!("structure_consistent_with:{inferred}")),
            Some(false) => {
                warnings.push(format!("expected_structure_missing_for:{inferred}"));
                confidence = Confidence::Low;
            }
            None => {}
        }
        return Disposition {
            kind: DispositionKind::KnownTe,
            label: inferred.clone(),
            confidence,
            evidence,
            warnings,
        };
    }

    // 3. Tandem repeats and low complexity
    if let Some(tandem) = structure
        .tandem
        .as_ref()
        .filter(|t2| t2.copies >= t.min_tandem_copies)
    {
        let (kind, label) = if tandem.period <= t.simple_repeat_max_period {
            (DispositionKind::SimpleRepeat, "Simple_repeat")
        } else {
            (DispositionKind::Satellite, "Satellite")
        };
        return Disposition {
            kind,
            label: label.into(),
            confidence: if tandem.identity >= 0.85 {
                Confidence::High
            } else {
                Confidence::Medium
            },
            evidence,
            warnings,
        };
    }
    if structure.composition.low_complexity {
        return Disposition {
            kind: DispositionKind::LowComplexity,
            label: "Low_complexity".into(),
            confidence: Confidence::Medium,
            evidence,
            warnings,
        };
    }

    // 4. TE hallmarks without homology: candidate novel elements
    let coding = structure.longest_orf_aa() >= t.coding_orf_aa;
    let structural_label = if structure.ltr.is_some() {
        Some(("LTR/Unknown", Confidence::Medium))
    } else if structure.tir.is_some() {
        Some(("DNA/Unknown", Confidence::Medium))
    } else if structure.tail.is_some() && coding {
        Some(("LINE/Unknown", Confidence::Low))
    } else if structure.helitron_termini {
        Some(("RC/Unknown", Confidence::Low))
    } else if coding {
        warnings.push("coding_without_te_homology_may_be_host_gene".into());
        Some(("Unknown", Confidence::Low))
    } else {
        None
    };
    if let Some((label, confidence)) = structural_label {
        return Disposition {
            kind: DispositionKind::NovelCandidate,
            label: label.into(),
            confidence,
            evidence,
            warnings,
        };
    }

    // 5. Nothing explains the family
    evidence.push(format!(
        "no_protein_homology;no_terminal_repeats;no_tail;longest_orf_aa={}",
        structure.longest_orf_aa()
    ));
    Disposition {
        kind: DispositionKind::Unresolved,
        label: "Unknown".into(),
        confidence: Confidence::Low,
        evidence,
        warnings,
    }
}

/// Refine a disposition with evidence from recovered genomic copies:
/// TSD length versus the class expectation, copy number, and flank conservation.
pub fn apply_copy_evidence(mut d: Disposition, c: &CopyEvidence, length: usize) -> Disposition {
    d.evidence.push(format!(
        "copies:loci={},full_length={},identity={:.1}",
        c.loci, c.full_length, c.mean_identity
    ));
    let conserved_flanks = c.flank_similarity.is_some_and(|f| f >= 0.6);
    if conserved_flanks {
        d.warnings.push(format!(
            "conserved_flanks:{:.2}(truncated_consensus_or_segmental_duplication)",
            c.flank_similarity.unwrap_or_default()
        ));
    }
    if c.loci < 3 {
        d.warnings.push(format!("low_copy_number:{}", c.loci));
    }
    if let Some(tsd) = &c.tsd {
        d.evidence.push(format!(
            "tsd:{}bp({}),support={:.2},null={:.2},n={}",
            tsd.length, tsd.motif, tsd.support, tsd.null_support, tsd.copies_evaluated
        ));
    }
    match d.kind {
        DispositionKind::KnownTe | DispositionKind::Chimeric => {
            if let (Some(tsd), Some(expected)) = (&c.tsd, expected_tsd(&d.label)) {
                if expected.contains(&tsd.length) {
                    d.evidence.push(format!("tsd_consistent_with:{}", d.label));
                } else if tsd.support >= 0.6 && tsd.copies_evaluated >= 5 {
                    d.warnings
                        .push(format!("tsd_inconsistent_with:{}", d.label));
                    d.confidence = d.confidence.lower();
                }
            }
        }
        DispositionKind::NovelCandidate => {
            if let Some(tsd) = &c.tsd {
                let hints = superfamilies_with_tsd(top_level(&d.label), tsd.length);
                if !hints.is_empty() {
                    d.evidence.push(format!("tsd_suggests:{}", hints.join(",")));
                }
            }
            if c.loci < 3 || conserved_flanks {
                d.confidence = d.confidence.lower();
            }
        }
        // tRNAs and other short multicopy genes are excluded by length
        DispositionKind::Unresolved
            if c.loci >= 10 && c.full_length >= 3 && !conserved_flanks && length >= 120 =>
        {
            // a real dispersed repeat; only a TE hallmark (here, a TSD) makes it a
            // novel TE candidate, since every RepeatModeler family is multicopy
            if c.tsd.is_some() {
                d.kind = DispositionKind::NovelCandidate;
                d.confidence = Confidence::Low;
                d.evidence
                    .push("dispersed_multicopy_repeat_with_tsd".into());
            } else {
                d.kind = DispositionKind::InterspersedRepeat;
                d.confidence = Confidence::Medium;
                d.evidence.push("dispersed_multicopy_repeat".into());
            }
        }
        _ => {}
    }
    d
}

/// Keep a specific non-TE label from the original classifier (tRNA, rRNA,
/// satellite, ...) unless TE proteins or domains contradict it.
pub fn respect_original_label(
    mut d: Disposition,
    original: Option<&str>,
    te_coding_evidence: bool,
) -> Disposition {
    let Some(original) = original.filter(|o| is_non_te(o)) else {
        return d;
    };
    let non_te_kind = matches!(
        d.kind,
        DispositionKind::HostOrArtefact
            | DispositionKind::Satellite
            | DispositionKind::SimpleRepeat
            | DispositionKind::LowComplexity
    );
    if non_te_kind || te_coding_evidence {
        return d;
    }
    d.warnings.push(format!(
        "kept_original_label:{original};conflicting:{}",
        d.label
    ));
    d.kind = match top_level(original) {
        "Satellite" => DispositionKind::Satellite,
        "Simple_repeat" => DispositionKind::SimpleRepeat,
        "Low_complexity" => DispositionKind::LowComplexity,
        _ => DispositionKind::HostOrArtefact,
    };
    d.label = original.to_owned();
    d.confidence = Confidence::Medium;
    d
}

fn structural_notes(s: &StructuralEvidence, t: &DecisionThresholds) -> Vec<String> {
    let mut notes = Vec::new();
    if let Some(r) = &s.ltr {
        notes.push(format!(
            "terminal_direct_repeat:{}bp@{:.2}",
            r.length, r.identity
        ));
    }
    if let Some(r) = &s.tir {
        notes.push(format!(
            "terminal_inverted_repeat:{}bp@{:.2}",
            r.length, r.identity
        ));
    }
    if let Some(tail) = s.tail {
        notes.push(tail.as_str().into());
    }
    if let Some(tandem) = &s.tandem {
        notes.push(format!(
            "tandem:period={}@{:.2}x{:.1}",
            tandem.period, tandem.identity, tandem.copies
        ));
    }
    if s.tg_ca_termini {
        notes.push("TG..CA_termini".into());
    }
    if s.helitron_termini {
        notes.push("TC..CTRR_termini".into());
    }
    if s.longest_orf_aa() >= t.coding_orf_aa {
        notes.push(format!("orf:{}aa", s.longest_orf_aa()));
    }
    notes
}

#[cfg(test)]
mod tests {
    use super::super::domains::{interpret, role_of, DomainHit};
    use super::super::structure::{analyse, StructureParams, TerminalRepeat};
    use super::*;

    fn random(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize]
            })
            .collect()
    }

    fn support(class: &str, bits: f64, cov: f64) -> HomologySummary {
        HomologySummary {
            classes: vec![ClassSupport {
                class: class.into(),
                best_bitscore: bits,
                best_subject: "s".into(),
                best_identity: 50.0,
                best_subject_coverage: 0.9,
                query_coverage: cov,
                n_hits: 1,
                intervals: vec![(0, 100)],
            }],
            ..Default::default()
        }
    }

    #[test]
    fn homology_with_consistent_structure_is_high_confidence() {
        let mut s = analyse(&random(4000, 1), &StructureParams::default());
        s.ltr = Some(TerminalRepeat {
            length: 300,
            identity: 0.95,
            five_prime_offset: 0,
            three_prime_offset: 0,
        });
        let d = decide(
            &support("LTR/Gypsy", 400.0, 0.4),
            &HomologySummary::default(),
            &DomainEvidence::default(),
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(
            (d.kind, d.label.as_str(), d.confidence),
            (DispositionKind::KnownTe, "LTR/Gypsy", Confidence::High)
        );
    }

    #[test]
    fn missing_expected_structure_lowers_confidence() {
        let s = analyse(&random(4000, 2), &StructureParams::default());
        let d = decide(
            &support("DNA/hAT", 400.0, 0.4),
            &HomologySummary::default(),
            &DomainEvidence::default(),
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(d.confidence, Confidence::Medium);
        assert!(d.warnings[0].starts_with("expected_structure_missing"));
    }

    #[test]
    fn satellites_and_simple_repeats() {
        let t = DecisionThresholds::default();
        let sat = analyse(&random(171, 3).repeat(8), &StructureParams::default());
        assert_eq!(
            decide(
                &HomologySummary::default(),
                &HomologySummary::default(),
                &DomainEvidence::default(),
                &sat,
                &t
            )
            .kind,
            DispositionKind::Satellite
        );
        let simple = analyse(&b"CAGT".repeat(50), &StructureParams::default());
        assert_eq!(
            decide(
                &HomologySummary::default(),
                &HomologySummary::default(),
                &DomainEvidence::default(),
                &simple,
                &t
            )
            .kind,
            DispositionKind::SimpleRepeat
        );
    }

    #[test]
    fn structural_hallmarks_without_homology_are_novel_candidates() {
        let t = DecisionThresholds::default();
        let ltr = random(250, 4);
        let seq = [ltr.clone(), random(3000, 5), ltr].concat();
        let d = decide(
            &HomologySummary::default(),
            &HomologySummary::default(),
            &DomainEvidence::default(),
            &analyse(&seq, &StructureParams::default()),
            &t,
        );
        assert_eq!(
            (d.kind, d.label.as_str()),
            (DispositionKind::NovelCandidate, "LTR/Unknown")
        );
        // a short non-coding poly-A family may be a SINE or a truncated LINE 3' end
        let sine = [random(300, 6), b"A".repeat(20)].concat();
        let d = decide(
            &HomologySummary::default(),
            &HomologySummary::default(),
            &DomainEvidence::default(),
            &analyse(&sine, &StructureParams::default()),
            &t,
        );
        assert_eq!(
            (d.kind, d.label.as_str()),
            (DispositionKind::Unresolved, "Unknown")
        );
        assert!(d.evidence.contains(&"polyA_3prime".to_string()));
    }

    #[test]
    fn unexplained_families_are_unresolved_with_reason() {
        let d = decide(
            &HomologySummary::default(),
            &HomologySummary::default(),
            &DomainEvidence::default(),
            &analyse(&random(2000, 7), &StructureParams::default()),
            &DecisionThresholds::default(),
        );
        assert_eq!(d.kind, DispositionKind::Unresolved);
        assert!(d.evidence.last().unwrap().contains("no_protein_homology"));
    }

    #[test]
    fn nucleotide_only_homology_classifies_non_coding_families() {
        let s = analyse(
            &[random(300, 8), b"A".repeat(20)].concat(),
            &StructureParams::default(),
        );
        let d = decide(
            &HomologySummary::default(),
            &support("SINE/tRNA", 400.0, 0.9),
            &DomainEvidence::default(),
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(
            (d.kind, d.label.as_str(), d.confidence),
            (DispositionKind::KnownTe, "SINE/tRNA", Confidence::High)
        );
        assert!(d.evidence[0].starts_with("nucleotide:SINE/tRNA"));
    }

    #[test]
    fn non_te_references_are_host_or_artefact() {
        let s = analyse(&random(2000, 9), &StructureParams::default());
        let d = decide(
            &HomologySummary::default(),
            &support("rRNA", 900.0, 0.95),
            &DomainEvidence::default(),
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(
            (d.kind, d.label.as_str()),
            (DispositionKind::HostOrArtefact, "rRNA")
        );
        // weak coverage is ignored
        let d = decide(
            &HomologySummary::default(),
            &support("rRNA", 900.0, 0.05),
            &DomainEvidence::default(),
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(d.kind, DispositionKind::Unresolved);
    }

    #[test]
    fn recurrences_of_our_novel_families_stay_novel() {
        let s = analyse(&random(2000, 10), &StructureParams::default());
        let d = decide(
            &HomologySummary::default(),
            &support("LTR/Unknown", 900.0, 0.9),
            &DomainEvidence::default(),
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(
            (d.kind, d.label.as_str()),
            (DispositionKind::NovelCandidate, "LTR/Unknown")
        );
        assert!(d
            .evidence
            .iter()
            .any(|e| e.starts_with("recurs_with_reference_family")));
    }

    #[test]
    fn protein_and_nucleotide_conflicts_lower_confidence() {
        let mut s = analyse(&random(4000, 11), &StructureParams::default());
        s.ltr = Some(TerminalRepeat {
            length: 300,
            identity: 0.95,
            five_prime_offset: 0,
            three_prime_offset: 0,
        });
        let d = decide(
            &support("LTR/Gypsy", 400.0, 0.4),
            &support("DNA/hAT", 500.0, 0.8),
            &DomainEvidence::default(),
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(
            (d.kind, d.confidence),
            (DispositionKind::KnownTe, Confidence::Medium)
        );
        assert!(d
            .warnings
            .iter()
            .any(|w| w.starts_with("nucleotide_class_conflict")));
    }

    fn domain_evidence(domains: &[(&str, usize)]) -> DomainEvidence {
        let hits = domains
            .iter()
            .map(|(d, start)| DomainHit {
                family_id: "f".into(),
                domain: (*d).into(),
                role: role_of(d).unwrap(),
                orf: "f|orf1|+:0-9000".into(),
                strand: '+',
                score: 50.0,
                i_evalue: 1e-20,
                coding_start: *start,
                coding_end: start + 300,
            })
            .collect::<Vec<_>>();
        interpret(&hits.iter().collect::<Vec<_>>())
    }

    #[test]
    fn domain_architecture_classifies_without_homology() {
        let s = analyse(&random(5000, 12), &StructureParams::default());
        let d = decide(
            &HomologySummary::default(),
            &HomologySummary::default(),
            &domain_evidence(&[("RVT_1", 0), ("RNase_H", 700), ("rve", 1400)]),
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(
            (d.kind, d.label.as_str()),
            (DispositionKind::KnownTe, "LTR/Gypsy")
        );
        assert!(d.evidence.contains(&"domains:RT-RH-INT".to_string()));
        // no LTRs in a random sequence
        assert_eq!(d.confidence, Confidence::Low);
    }

    #[test]
    fn novel_domain_architectures_are_candidates_even_with_partial_homology() {
        let s = analyse(&random(5000, 13), &StructureParams::default());
        let novel = domain_evidence(&[("RVT_1", 0), ("Dimer_Tnp_hAT", 1200)]);
        let d = decide(
            &HomologySummary::default(),
            &HomologySummary::default(),
            &novel,
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(
            (d.kind, d.label.as_str()),
            (DispositionKind::NovelCandidate, "Unknown")
        );
        assert!(d.evidence[0].starts_with("novel_domain_architecture:RT+DDE-transposase"));
        // with homology the class stands, but the oddity is flagged
        let d = decide(
            &support("LINE/L1", 400.0, 0.4),
            &HomologySummary::default(),
            &novel,
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(d.kind, DispositionKind::KnownTe);
        assert!(d
            .warnings
            .iter()
            .any(|w| w.starts_with("novel_domain_architecture")));
    }

    #[test]
    fn domain_conflicts_lower_homology_confidence() {
        let mut s = analyse(&random(5000, 14), &StructureParams::default());
        s.ltr = Some(TerminalRepeat {
            length: 300,
            identity: 0.95,
            five_prime_offset: 0,
            three_prime_offset: 0,
        });
        let copia = domain_evidence(&[("rve", 0), ("RVT_2", 900)]);
        let d = decide(
            &support("LTR/Gypsy", 400.0, 0.4),
            &HomologySummary::default(),
            &copia,
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(
            d.confidence,
            Confidence::High,
            "same top-level class: no conflict"
        );
        let dna = domain_evidence(&[("Dimer_Tnp_hAT", 0)]);
        let d = decide(
            &support("LTR/Gypsy", 400.0, 0.4),
            &HomologySummary::default(),
            &dna,
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(d.confidence, Confidence::Medium);
        assert!(d
            .warnings
            .iter()
            .any(|w| w == "domain_class_conflict:DNA/hAT"));
    }

    fn copies(loci: usize, full: usize, tsd: Option<(usize, f64)>, flanks: f64) -> CopyEvidence {
        CopyEvidence {
            loci,
            full_length: full,
            mean_identity: 90.0,
            tsd: tsd.map(|(length, support)| super::super::copies::Tsd {
                length,
                support,
                null_support: 0.1,
                motif: "TA".into(),
                copies_evaluated: 10,
            }),
            flank_similarity: Some(flanks),
            genome_bp: 0,
        }
    }

    fn disposition(kind: DispositionKind, label: &str, confidence: Confidence) -> Disposition {
        Disposition {
            kind,
            label: label.into(),
            confidence,
            evidence: vec![],
            warnings: vec![],
        }
    }

    #[test]
    fn tsd_confirms_or_contradicts_known_classes() {
        let d = apply_copy_evidence(
            disposition(
                DispositionKind::KnownTe,
                "DNA/TcMar-Mariner",
                Confidence::High,
            ),
            &copies(20, 8, Some((2, 0.8)), 0.3),
            500,
        );
        assert!(d
            .evidence
            .iter()
            .any(|e| e == "tsd_consistent_with:DNA/TcMar-Mariner"));
        assert_eq!(d.confidence, Confidence::High);
        let d = apply_copy_evidence(
            disposition(DispositionKind::KnownTe, "DNA/hAT", Confidence::High),
            &copies(20, 8, Some((2, 0.8)), 0.3),
            500,
        );
        assert_eq!(d.confidence, Confidence::Medium);
        assert!(d
            .warnings
            .iter()
            .any(|w| w == "tsd_inconsistent_with:DNA/hAT"));
    }

    #[test]
    fn tsd_hints_superfamilies_for_novel_candidates() {
        let d = apply_copy_evidence(
            disposition(
                DispositionKind::NovelCandidate,
                "DNA/Unknown",
                Confidence::Medium,
            ),
            &copies(30, 10, Some((8, 0.7)), 0.3),
            500,
        );
        assert!(d
            .evidence
            .iter()
            .any(|e| e == "tsd_suggests:DNA/hAT,DNA/MULE,DNA/P,DNA/Merlin"));
        assert_eq!(d.confidence, Confidence::Medium);
        let d = apply_copy_evidence(
            disposition(
                DispositionKind::NovelCandidate,
                "DNA/Unknown",
                Confidence::Medium,
            ),
            &copies(2, 1, None, 0.9),
            500,
        );
        assert_eq!(d.confidence, Confidence::Low);
        assert!(d.warnings.iter().any(|w| w.starts_with("conserved_flanks")));
        assert!(d.warnings.iter().any(|w| w == "low_copy_number:2"));
    }

    #[test]
    fn dispersed_multicopy_unresolved_families_become_candidates() {
        let d = apply_copy_evidence(
            disposition(DispositionKind::Unresolved, "Unknown", Confidence::Low),
            &copies(40, 6, None, 0.3),
            500,
        );
        assert_eq!(d.kind, DispositionKind::InterspersedRepeat);
        let d = apply_copy_evidence(
            disposition(DispositionKind::Unresolved, "Unknown", Confidence::Low),
            &copies(40, 6, Some((5, 0.7)), 0.3),
            500,
        );
        assert_eq!(d.kind, DispositionKind::NovelCandidate);
        let d = apply_copy_evidence(
            disposition(DispositionKind::Unresolved, "Unknown", Confidence::Low),
            &copies(40, 6, None, 0.8),
            500,
        );
        assert_eq!(
            d.kind,
            DispositionKind::Unresolved,
            "segmental-duplication-like flanks"
        );
    }

    #[test]
    fn short_dispersed_families_are_not_promoted() {
        let d = apply_copy_evidence(
            disposition(DispositionKind::Unresolved, "Unknown", Confidence::Low),
            &copies(40, 30, None, 0.3),
            90,
        );
        assert_eq!(d.kind, DispositionKind::Unresolved);
    }

    #[test]
    fn original_non_te_labels_are_kept_without_te_coding_evidence() {
        let novel = disposition(DispositionKind::NovelCandidate, "Unknown", Confidence::Low);
        let d = respect_original_label(novel.clone(), Some("tRNA"), false);
        assert_eq!(
            (d.kind, d.label.as_str()),
            (DispositionKind::HostOrArtefact, "tRNA")
        );
        assert!(d.warnings[0].starts_with("kept_original_label:tRNA"));
        assert_eq!(
            respect_original_label(novel.clone(), Some("tRNA"), true).kind,
            DispositionKind::NovelCandidate
        );
        assert_eq!(
            respect_original_label(novel, Some("LTR/Gypsy"), false).kind,
            DispositionKind::NovelCandidate
        );
    }

    #[test]
    fn short_family_matching_a_sliver_of_a_long_reference_is_low_confidence() {
        let s = analyse(&random(132, 15), &StructureParams::default());
        let mut nucl = support("LINE/L1", 124.0, 0.77);
        nucl.classes[0].best_subject_coverage = 0.02;
        let d = decide(
            &HomologySummary::default(),
            &nucl,
            &DomainEvidence::default(),
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(d.confidence, Confidence::Low);
        assert!(d
            .warnings
            .iter()
            .any(|w| w.starts_with("short_match_within_long_reference")));
    }

    #[test]
    fn weak_nucleotide_homology_is_a_low_confidence_call() {
        let s = analyse(&random(1500, 16), &StructureParams::default());
        let d = decide(
            &HomologySummary::default(),
            &support("LINE/RTE", 80.0, 0.2),
            &DomainEvidence::default(),
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(
            (d.kind, d.label.as_str(), d.confidence),
            (DispositionKind::KnownTe, "LINE/RTE", Confidence::Low)
        );
        assert!(d.warnings.contains(&"weak_nucleotide_homology".to_string()));
    }

    #[test]
    fn composite_sines_and_helitrons_take_their_own_class() {
        let s = analyse(&random(600, 17), &StructureParams::default());
        let mut nucl = support("LINE/L2", 300.0, 0.5);
        nucl.classes.push(ClassSupport {
            class: "SINE/tRNA".into(),
            best_bitscore: 200.0,
            best_subject: "x".into(),
            best_identity: 80.0,
            best_subject_coverage: 0.9,
            query_coverage: 0.4,
            n_hits: 1,
            intervals: vec![(0, 200)],
        });
        nucl.chimeric_classes = vec!["LINE".into(), "SINE".into()];
        let d = decide(
            &HomologySummary::default(),
            &nucl,
            &DomainEvidence::default(),
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(
            (d.kind, d.label.as_str(), d.confidence),
            (DispositionKind::KnownTe, "SINE/tRNA", Confidence::Low)
        );
        assert!(d
            .warnings
            .iter()
            .any(|w| w == "composite_element:SINE/tRNA"));
        nucl.classes[1].class = "DNA/hAT".into();
        nucl.chimeric_classes = vec!["LINE".into(), "DNA".into()];
        let d = decide(
            &HomologySummary::default(),
            &nucl,
            &DomainEvidence::default(),
            &s,
            &DecisionThresholds::default(),
        );
        assert_eq!(d.kind, DispositionKind::Chimeric);
    }
}
