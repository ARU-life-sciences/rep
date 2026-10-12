//! Combine homology and structural evidence into a per-family disposition.
//!
//! The rules are deliberately transparent: every disposition lists the
//! evidence that produced it, and anything not explained is `unresolved` with
//! the reason recorded.

use serde::Serialize;

use super::{
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
            nucl_min_bitscore: 100.0,
            nucl_min_coverage: 0.3,
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
    structure: &StructuralEvidence,
    t: &DecisionThresholds,
) -> Disposition {
    let mut evidence = structural_notes(structure, t);
    let mut warnings = Vec::new();
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

    // 2. Tandem repeats and low complexity
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

    // 3. TE hallmarks without homology: candidate novel elements
    let coding = structure.longest_orf_aa() >= t.coding_orf_aa;
    let structural_label = if structure.ltr.is_some() {
        Some(("LTR/Unknown", Confidence::Medium))
    } else if structure.tir.is_some() {
        Some(("DNA/Unknown", Confidence::Medium))
    } else if structure.tail.is_some() && coding {
        Some(("LINE/Unknown", Confidence::Low))
    } else if structure.tail.is_some() && structure.length <= t.sine_max_length {
        Some(("SINE/Unknown", Confidence::Low))
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

    // 4. Nothing explains the family
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
            &analyse(&seq, &StructureParams::default()),
            &t,
        );
        assert_eq!(
            (d.kind, d.label.as_str()),
            (DispositionKind::NovelCandidate, "LTR/Unknown")
        );
        let sine = [random(300, 6), b"A".repeat(20)].concat();
        let d = decide(
            &HomologySummary::default(),
            &HomologySummary::default(),
            &analyse(&sine, &StructureParams::default()),
            &t,
        );
        assert_eq!(d.label, "SINE/Unknown");
    }

    #[test]
    fn unexplained_families_are_unresolved_with_reason() {
        let d = decide(
            &HomologySummary::default(),
            &HomologySummary::default(),
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
            &support("rRNA", 900.0, 0.1),
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
}
