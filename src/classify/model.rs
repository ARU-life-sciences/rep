//! The learned top-level classifier (ROADMAP M8c).
//!
//! A multinomial logistic regression over the `features.tsv` values, fitted
//! offline by `scripts/classifier/fit.py` and shipped as `build/classifier.json`
//! in the refdb. It runs beside the rule cascade: the cascade's non-TE calls
//! (simple repeats, satellites, RNA and host genes) take precedence, otherwise
//! the model picks the class and the superfamily comes from the best homology
//! hit within that class.
//!
//! Which call is final (M8d, `scripts/classifier/switch_report.py`): the model's
//! at [`LOW`] or above where the cascade made a TE call, and at [`MEDIUM`] or above
//! where it made none (so Unknowns are only labelled at medium or high), except that
//! the cascade keeps its non-TE calls, its calls backed by a protein hit in their
//! own class, and families the model calls Retroposon or nonTE (held out, the model
//! was worse than the cascade there).

use std::fs;
use std::path::Path;

use serde::{Deserialize, Serialize};

use super::{
    decide::{is_non_te, Confidence, Disposition, DispositionKind},
    domains::DomainEvidence,
    features::{group, numeric_header, FEATURES_SCHEMA},
    homology::HomologySummary,
};
use crate::{Error, ErrorKind, Result};

pub const CLASSIFIER_SCHEMA: &str = "rep.classifier.v1";

/// Posterior bands, from held-out agreement with RepeatModeler labels:
/// ≥0.99 agreed 98.2%, 0.90–0.99 about 96%, 0.70–0.90 about 85%.
const HIGH: f64 = 0.99;
const MEDIUM: f64 = 0.9;
const LOW: f64 = 0.7;
/// Model classes that never replace the cascade's call.
const KEEP_CASCADE: &[&str] = &["Retroposon", "nonTE"];

#[derive(Deserialize)]
struct TransformFile {
    columns: Vec<String>,
    log1p: Vec<String>,
    missing_indicators: Vec<String>,
    mean: Vec<f64>,
    std: Vec<f64>,
}

#[derive(Deserialize)]
struct ModelFile {
    schema: String,
    features_schema: String,
    classes: Vec<String>,
    transform: TransformFile,
    weights: Vec<Vec<f64>>,
    bias: Vec<f64>,
}

#[derive(Debug)]
pub struct Classifier {
    classes: Vec<String>,
    /// Position in [`numeric_header`] of each model column.
    index: Vec<usize>,
    log: Vec<bool>,
    /// Model columns that also get a missing-value indicator.
    missing: Vec<usize>,
    mean: Vec<f64>,
    std: Vec<f64>,
    weights: Vec<Vec<f64>>,
    bias: Vec<f64>,
}

/// The model's call for one family.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ModelCall {
    /// Top-level class: one of the model's classes, e.g. `DNA` or `nonTE`.
    pub class: String,
    pub label: String,
    pub probability: f64,
    /// `high`, `medium`, `low`, or `none` below the lowest band.
    pub confidence: String,
    /// The cascade's non-TE call was kept instead of the model's.
    pub gated: bool,
    /// The model's call replaces the cascade's (see the module docs).
    #[serde(default)]
    pub decides: bool,
}

fn invalid(msg: String) -> Error {
    Error::new(ErrorKind::GenericCli(format!("classifier: {msg}")))
}

impl Classifier {
    pub fn load(path: &Path) -> Result<Self> {
        let text = fs::read_to_string(path)?;
        Self::from_json(&text).map_err(|e| invalid(format!("{}: {e}", path.display())))
    }

    pub fn from_json(text: &str) -> std::result::Result<Self, String> {
        let m: ModelFile = serde_json::from_str(text).map_err(|e| e.to_string())?;
        if m.schema != CLASSIFIER_SCHEMA {
            return Err(format!(
                "schema {} (expected {CLASSIFIER_SCHEMA})",
                m.schema
            ));
        }
        if m.features_schema != FEATURES_SCHEMA {
            return Err(format!(
                "trained on features {} (this rep writes {FEATURES_SCHEMA})",
                m.features_schema
            ));
        }
        let header = numeric_header();
        let t = &m.transform;
        let index = t
            .columns
            .iter()
            .map(|c| {
                header
                    .iter()
                    .position(|h| h == c)
                    .ok_or_else(|| format!("unknown feature column {c}"))
            })
            .collect::<std::result::Result<Vec<_>, _>>()?;
        let missing = t
            .missing_indicators
            .iter()
            .map(|c| {
                t.columns
                    .iter()
                    .position(|h| h == c)
                    .ok_or_else(|| format!("missing indicator for unknown column {c}"))
            })
            .collect::<std::result::Result<Vec<_>, _>>()?;
        let d = t.columns.len() + missing.len();
        if t.mean.len() != d || t.std.len() != d {
            return Err(format!(
                "transform has {} means for {d} features",
                t.mean.len()
            ));
        }
        if m.weights.len() != m.classes.len()
            || m.bias.len() != m.classes.len()
            || m.weights.iter().any(|w| w.len() != d)
        {
            return Err("weights do not match the classes and features".into());
        }
        Ok(Self {
            log: t.columns.iter().map(|c| t.log1p.contains(c)).collect(),
            classes: m.classes,
            index,
            missing,
            mean: t.mean.clone(),
            std: t.std.clone(),
            weights: m.weights,
            bias: m.bias,
        })
    }

    pub fn classes(&self) -> &[String] {
        &self.classes
    }

    /// Class probabilities for one family's [`super::features::values`].
    pub fn probabilities(&self, values: &[Option<f64>]) -> Vec<f64> {
        let mut x = self
            .index
            .iter()
            .zip(&self.log)
            .map(|(&i, &log)| {
                let v = values[i].unwrap_or(0.0);
                if log {
                    v.max(0.0).ln_1p()
                } else {
                    v
                }
            })
            .collect::<Vec<_>>();
        x.extend(self.missing.iter().map(|&j| {
            if values[self.index[j]].is_none() {
                1.0
            } else {
                0.0
            }
        }));
        for ((v, m), s) in x.iter_mut().zip(&self.mean).zip(&self.std) {
            *v = (*v - m) / s;
        }
        let logits = self
            .weights
            .iter()
            .zip(&self.bias)
            .map(|(w, b)| b + w.iter().zip(&x).map(|(a, b)| a * b).sum::<f64>())
            .collect::<Vec<_>>();
        let max = logits.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let exp = logits.iter().map(|z| (z - max).exp()).collect::<Vec<_>>();
        let total: f64 = exp.iter().sum();
        exp.into_iter().map(|e| e / total).collect()
    }

    /// The model's call, with the cascade's non-TE calls kept.
    pub fn call(
        &self,
        values: &[Option<f64>],
        cascade: &Disposition,
        protein: &HomologySummary,
        nucleotide: &HomologySummary,
        domains: &DomainEvidence,
    ) -> ModelCall {
        if cascade_non_te(cascade) {
            return ModelCall {
                class: "nonTE".into(),
                label: cascade.label.clone(),
                probability: 1.0,
                confidence: cascade.confidence.as_str().into(),
                gated: true,
                decides: false,
            };
        }
        let p = self.probabilities(values);
        let (k, &probability) = p
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .expect("a classifier has classes");
        let class = &self.classes[k];
        ModelCall {
            class: class.clone(),
            label: label_within(class, cascade, protein, nucleotide, domains),
            probability,
            confidence: band(probability).into(),
            gated: false,
            decides: probability
                >= if cascade_te_call(cascade) {
                    LOW
                } else {
                    MEDIUM
                }
                && !KEEP_CASCADE.contains(&class.as_str())
                && !protein_backed(values, cascade),
        }
    }
}

/// The model's call as the family's disposition: its label and band, with the
/// cascade's call kept as evidence. The disposition kind stays the cascade's.
pub fn apply(mut d: Disposition, m: &ModelCall) -> Disposition {
    d.evidence.push(format!(
        "model:{}@{:.3};cascade:{}",
        m.label, m.probability, d.label
    ));
    d.label = m.label.clone();
    d.confidence = match m.confidence.as_str() {
        "high" => Confidence::High,
        "medium" => Confidence::Medium,
        _ => Confidence::Low,
    };
    d
}

/// The cascade called a class (non-TE calls are gated before this matters).
fn cascade_te_call(d: &Disposition) -> bool {
    matches!(
        d.kind,
        DispositionKind::KnownTe | DispositionKind::Chimeric | DispositionKind::HostOrArtefact
    )
}

/// A cascade TE call with a protein (blastx) hit in its own class.
fn protein_backed(values: &[Option<f64>], cascade: &Disposition) -> bool {
    if !cascade_te_call(cascade) {
        return false;
    }
    let Some(class) = group(&cascade.label) else {
        return false;
    };
    let column = format!("prot_bits_{class}");
    numeric_header()
        .iter()
        .position(|h| *h == column)
        .and_then(|i| values[i])
        .is_some_and(|bits| bits > 0.0)
}

fn band(p: f64) -> &'static str {
    if p >= HIGH {
        "high"
    } else if p >= MEDIUM {
        "medium"
    } else if p >= LOW {
        "low"
    } else {
        "none"
    }
}

/// Simple repeats, satellites and non-TE references called by the cascade.
fn cascade_non_te(d: &Disposition) -> bool {
    match d.kind {
        DispositionKind::SimpleRepeat
        | DispositionKind::Satellite
        | DispositionKind::LowComplexity => true,
        DispositionKind::KnownTe | DispositionKind::Chimeric | DispositionKind::HostOrArtefact => {
            is_non_te(&d.label)
        }
        _ => false,
    }
}

/// The most specific label within `class`: the cascade's when it agrees, else
/// the best protein hit, the best nucleotide hit, or the domain-implied class.
fn label_within(
    class: &str,
    cascade: &Disposition,
    protein: &HomologySummary,
    nucleotide: &HomologySummary,
    domains: &DomainEvidence,
) -> String {
    let in_class = |label: &str| group(label) == Some(class);
    if in_class(&cascade.label) && cascade.label.contains('/') {
        return cascade.label.clone();
    }
    for summary in [protein, nucleotide] {
        // classes are in descending best-bitscore order
        if let Some(best) = summary.classes.iter().find(|c| in_class(&c.class)) {
            return best.class.clone();
        }
    }
    match domains.inferred.as_deref() {
        Some(d) if in_class(d) => d.to_owned(),
        _ if in_class(&cascade.label) => cascade.label.clone(),
        _ if class == "nonTE" => "Other".into(),
        _ => class.to_owned(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::classify::homology::ClassSupport;

    /// A two-class model on two features: protein DNA bits (log1p) and TSD length (optional).
    fn tiny() -> Classifier {
        let json = serde_json::json!({
            "schema": CLASSIFIER_SCHEMA,
            "features_schema": FEATURES_SCHEMA,
            "model": "multinomial_logistic",
            "classes": ["DNA", "LINE"],
            "transform": {
                "columns": ["prot_bits_DNA", "tsd_length"],
                "log1p": ["prot_bits_DNA"],
                "missing_indicators": ["tsd_length"],
                "features": ["prot_bits_DNA", "tsd_length", "tsd_length_missing"],
                "mean": [1.0, 0.0, 0.5],
                "std": [2.0, 1.0, 0.5]
            },
            "weights": [[2.0, 0.5, 0.0], [-2.0, 0.0, 1.0]],
            "bias": [0.0, 0.0]
        });
        Classifier::from_json(&json.to_string()).unwrap()
    }

    fn values(bits: f64, tsd: Option<f64>) -> Vec<Option<f64>> {
        let header = numeric_header();
        let mut v = vec![Some(0.0); header.len()];
        v[header.iter().position(|h| h == "prot_bits_DNA").unwrap()] = Some(bits);
        v[header.iter().position(|h| h == "tsd_length").unwrap()] = tsd;
        v
    }

    fn disposition(kind: DispositionKind, label: &str) -> Disposition {
        Disposition {
            kind,
            label: label.into(),
            confidence: Confidence::Medium,
            evidence: vec![],
            warnings: vec![],
        }
    }

    fn support(class: &str, bits: f64) -> ClassSupport {
        ClassSupport {
            class: class.into(),
            best_bitscore: bits,
            best_subject: "x".into(),
            best_identity: 80.0,
            best_subject_coverage: 0.5,
            query_coverage: 0.5,
            n_hits: 1,
            intervals: vec![],
        }
    }

    #[test]
    fn probabilities_follow_the_transform() {
        let m = tiny();
        // x = ((ln(1+e^5-1) - 1)/2, (8-0)/1, (0-0.5)/0.5) = (2, 8, -1)
        let p = m.probabilities(&values(5f64.exp() - 1.0, Some(8.0)));
        let (a, b) = (2.0 * 2.0 + 0.5 * 8.0, -2.0 * 2.0 - 1.0);
        let expect = 1.0 / (1.0 + (b - a as f64).exp());
        assert!((p[0] - expect).abs() < 1e-12, "{p:?}");
        // a missing TSD switches its indicator on: x = (2, 0, 1)
        let q = m.probabilities(&values(5f64.exp() - 1.0, None));
        let expect = 1.0 / (1.0 + ((-4.0 + 1.0) - 4.0f64).exp());
        assert!((q[0] - expect).abs() < 1e-12, "{q:?}");
    }

    #[test]
    fn model_class_takes_the_best_label_within_it() {
        let m = tiny();
        let protein = HomologySummary {
            classes: vec![support("LINE/RTE", 300.0), support("DNA/hAT", 120.0)],
            ..Default::default()
        };
        let cascade = disposition(DispositionKind::KnownTe, "LINE/RTE");
        let c = m.call(
            &values(400.0, Some(8.0)),
            &cascade,
            &protein,
            &HomologySummary::default(),
            &DomainEvidence::default(),
        );
        assert_eq!((c.class.as_str(), c.label.as_str()), ("DNA", "DNA/hAT"));
        assert!(!c.gated);
        // when the cascade's label is in the model's class it is kept
        let cascade = disposition(DispositionKind::KnownTe, "DNA/TcMar-Tc1");
        let c = m.call(
            &values(400.0, Some(8.0)),
            &cascade,
            &protein,
            &HomologySummary::default(),
            &DomainEvidence::default(),
        );
        assert_eq!(c.label, "DNA/TcMar-Tc1");
        assert_eq!(c.confidence, "high");
    }

    #[test]
    fn protein_backed_cascade_calls_and_kept_classes_are_not_replaced() {
        let m = tiny();
        let nothing = HomologySummary::default();
        let domains = DomainEvidence::default();
        // a LINE call without protein support: the confident DNA call decides
        let line = disposition(DispositionKind::KnownTe, "LINE/RTE");
        let c = m.call(
            &values(400.0, Some(8.0)),
            &line,
            &nothing,
            &nothing,
            &domains,
        );
        assert!(c.decides);
        let d = apply(line.clone(), &c);
        assert_eq!((d.label.as_str(), d.confidence), ("DNA", Confidence::High));
        assert_eq!(d.kind, DispositionKind::KnownTe);
        assert!(d.evidence.iter().any(|e| e.ends_with("cascade:LINE/RTE")));
        // a DNA call backed by a DNA protein hit stands
        let dna = disposition(DispositionKind::KnownTe, "DNA/hAT");
        assert!(
            !m.call(
                &values(400.0, Some(8.0)),
                &dna,
                &nothing,
                &nothing,
                &domains
            )
            .decides
        );
        // an unresolved family takes a confident call, but not an unconfident one
        let none = disposition(DispositionKind::Unresolved, "Unknown");
        assert!(
            m.call(
                &values(400.0, Some(8.0)),
                &none,
                &nothing,
                &nothing,
                &domains
            )
            .decides
        );
        // log1p(bits) = 1.5 puts the two classes level
        let c = m.call(
            &values(1.5f64.exp() - 1.0, None),
            &none,
            &nothing,
            &nothing,
            &domains,
        );
        assert!(c.probability < LOW && !c.decides, "{c:?}");
        // at p = 0.8 the model fixes an unsupported cascade call but leaves an Unknown alone
        // (logit gap 4z - 1 = ln 4 at standardised bits z, log1p(bits) = 1 + 2z)
        let z = (4f64.ln() + 1.0) / 4.0;
        let low = values((1.0 + 2.0 * z).exp() - 1.0, None);
        let c = m.call(&low, &none, &nothing, &nothing, &domains);
        assert!((c.probability - 0.8).abs() < 1e-9 && !c.decides, "{c:?}");
        assert!(m.call(&low, &line, &nothing, &nothing, &domains).decides);
        // the model's Retroposon and nonTE calls never decide
        let mut kept = tiny();
        kept.classes = vec!["Retroposon".into(), "nonTE".into()];
        let c = kept.call(
            &values(400.0, Some(8.0)),
            &none,
            &nothing,
            &nothing,
            &domains,
        );
        assert!(c.probability >= LOW && !c.decides);
    }

    #[test]
    fn cascade_non_te_calls_are_kept() {
        let m = tiny();
        for d in [
            disposition(DispositionKind::Satellite, "Satellite"),
            disposition(DispositionKind::HostOrArtefact, "rRNA"),
        ] {
            let c = m.call(
                &values(400.0, Some(8.0)),
                &d,
                &HomologySummary::default(),
                &HomologySummary::default(),
                &DomainEvidence::default(),
            );
            assert!(c.gated && !c.decides);
            assert_eq!(
                (c.class.as_str(), c.label.as_str()),
                ("nonTE", d.label.as_str())
            );
        }
    }

    #[test]
    fn mismatched_models_are_rejected() {
        let old = serde_json::json!({
            "schema": CLASSIFIER_SCHEMA, "features_schema": "rep.features.v0", "classes": ["DNA"],
            "transform": {"columns": [], "log1p": [], "missing_indicators": [], "mean": [], "std": []},
            "weights": [[]], "bias": [0.0]
        });
        assert!(Classifier::from_json(&old.to_string())
            .unwrap_err()
            .contains("rep.features.v0"));
        let unknown = serde_json::json!({
            "schema": CLASSIFIER_SCHEMA, "features_schema": FEATURES_SCHEMA, "classes": ["DNA"],
            "transform": {"columns": ["no_such"], "log1p": [], "missing_indicators": [], "mean": [0.0], "std": [1.0]},
            "weights": [[1.0]], "bias": [0.0]
        });
        assert!(Classifier::from_json(&unknown.to_string())
            .unwrap_err()
            .contains("no_such"));
    }
}
