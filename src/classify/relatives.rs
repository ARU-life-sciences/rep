//! Second pass for Unknowns: nucleotide homology to families already called
//! as TEs, in the same library or in related genomes (`--relatives`).
//!
//! Most Unknowns that nothing else explains match no reference element, but
//! many are pieces of an element that RepeatModeler split into several
//! consensi, or relatives of a family called in a related genome. Held out on
//! labelled families (40 Ensembl libraries, 2026-10-14), the class of the best
//! such match agreed with the family's own call 0.93–0.99 of the time where
//! it covered at least 60% of the family, and 0.81–0.95 below that, where the
//! confusions are the expected LINE/SINE (shared 3' ends) and RC/DNA ones.
//!
//! On the 17-moss panel (2026-10-19), with the protein-coding regions of
//! protein-labelled families removed so the queries resemble non-coding
//! Unknowns, the call was right 97% of the time, also below 20% coverage. It
//! was right only 67% of the time when another class scored nearly as well
//! (best score < 1.5 x the best other class's), which now makes no call.
//!
//! Only families called by the cascade or the model serve as references,
//! never families labelled by this pass or from genomic neighbours, so labels
//! do not travel along chains of inference.

use serde::Serialize;

use super::{
    decide::{is_non_te, Confidence, Disposition},
    homology::HomologySummary,
    taxonomy::top_level,
};

/// The best class's score over the best other class's, below which the
/// match does not tell the two apart.
pub const MIN_MARGIN: f64 = 1.5;

/// Prefix of references from the library being classified, to tell them
/// apart from `--relatives` references in the evidence.
pub const LIBRARY_PREFIX: &str = "library|";

/// A label from a related family.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct RelativeCall {
    pub label: String,
    pub subject: String,
    pub bitscore: f64,
    /// Fraction of the family covered by hits to the label's class.
    pub coverage: f64,
    pub confidence: Confidence,
}

/// The best related family's label, unless the match is too weak, two
/// classes compete for different parts of the family, or another class
/// scores nearly as well.
pub fn call(summary: &HomologySummary, medium_coverage: f64) -> Option<RelativeCall> {
    if !summary.chimeric_classes.is_empty() {
        return None;
    }
    let best = summary.best().filter(|b| !is_non_te(&b.class))?;
    let other = summary
        .classes
        .iter()
        .filter(|c| top_level(&c.class) != top_level(&best.class))
        .map(|c| c.best_bitscore)
        .fold(0.0, f64::max);
    if best.best_bitscore < MIN_MARGIN * other {
        return None;
    }
    Some(RelativeCall {
        label: best.class.clone(),
        subject: best.best_subject.clone(),
        bitscore: best.best_bitscore,
        coverage: best.query_coverage,
        confidence: if best.query_coverage >= medium_coverage {
            Confidence::Medium
        } else {
            Confidence::Low
        },
    })
}

/// The related family's label as the family's disposition; the kind stays.
pub fn apply(mut d: Disposition, r: &RelativeCall) -> Disposition {
    d.evidence.push(format!(
        "relative:{}@{:.0},{:.2};was:{}",
        r.subject, r.bitscore, r.coverage, d.label
    ));
    d.label = r.label.clone();
    d.confidence = r.confidence;
    d
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::classify::decide::DispositionKind;
    use crate::classify::homology::ClassSupport;

    fn support(class: &str, bits: f64, coverage: f64) -> ClassSupport {
        ClassSupport {
            class: class.into(),
            best_bitscore: bits,
            best_subject: format!("library|x#{class}"),
            best_identity: 80.0,
            best_subject_coverage: 0.5,
            query_coverage: coverage,
            n_hits: 1,
            intervals: vec![],
        }
    }

    #[test]
    fn coverage_sets_the_confidence() {
        let mut s = HomologySummary {
            classes: vec![
                support("LINE/L2", 300.0, 0.8),
                support("DNA/hAT", 90.0, 0.2),
            ],
            ..Default::default()
        };
        let r = call(&s, 0.6).unwrap();
        assert_eq!(
            (r.label.as_str(), r.confidence),
            ("LINE/L2", Confidence::Medium)
        );
        s.classes[0].query_coverage = 0.4;
        assert_eq!(call(&s, 0.6).unwrap().confidence, Confidence::Low);
        let d = apply(
            Disposition {
                kind: DispositionKind::InterspersedRepeat,
                label: "Unknown".into(),
                confidence: Confidence::Low,
                evidence: vec![],
                warnings: vec![],
            },
            &r,
        );
        assert_eq!(d.kind, DispositionKind::InterspersedRepeat);
        assert_eq!(d.label, "LINE/L2");
        assert_eq!(
            d.evidence,
            vec!["relative:library|x#LINE/L2@300,0.80;was:Unknown"]
        );
    }

    #[test]
    fn chimeric_and_non_te_matches_make_no_call() {
        let s = HomologySummary {
            classes: vec![
                support("LINE/L2", 300.0, 0.5),
                support("DNA/hAT", 250.0, 0.5),
            ],
            chimeric_classes: vec!["LINE".into(), "DNA".into()],
            ..Default::default()
        };
        assert!(call(&s, 0.6).is_none());
        let s = HomologySummary {
            classes: vec![support("Satellite", 300.0, 0.9)],
            ..Default::default()
        };
        assert!(call(&s, 0.6).is_none());
        assert!(call(&HomologySummary::default(), 0.6).is_none());
    }

    #[test]
    fn a_close_second_class_makes_no_call() {
        // overlapping hits, so not chimeric, but LTR nearly as good as LINE
        let mut s = HomologySummary {
            classes: vec![
                support("LINE/L1", 200.0, 0.5),
                support("LTR/Gypsy", 150.0, 0.5),
                support("LINE/RTE", 180.0, 0.5),
            ],
            ..Default::default()
        };
        assert!(call(&s, 0.6).is_none());
        // superfamilies of the same class do not compete
        s.classes[1].best_bitscore = 120.0;
        assert_eq!(call(&s, 0.6).unwrap().label, "LINE/L1");
    }
}
