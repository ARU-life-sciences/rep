//! Per-family feature table for a learned classifier (ROADMAP M8a).
//!
//! One row per classified family with every piece of evidence as a number,
//! so a model can weigh evidence instead of the rule cascade in `decide`.
//! Homology is summarised per top-level class and per source (BLASTX
//! proteins, superfamily profile HMMs, BLASTN consensi), so the model learns
//! how much each source is worth. Missing values (no copies, no TSD) are
//! empty fields; presence is also given as 0/1 columns where it matters.
//! The first columns identify the family and carry the original
//! (RepeatClassifier) label and our current call, as training targets and
//! for comparison; they are not features.

use std::collections::BTreeSet;
use std::fmt::Write as _;

use super::{
    copies::CopyEvidence,
    decide::is_non_te,
    domains::{DomainEvidence, Role},
    homology::HomologySummary,
    structure::{StructuralEvidence, Tail},
    taxonomy::canonical_class,
};

pub const FEATURES_SCHEMA: &str = "rep.features.v1";

/// Top-level classes scored separately; `nonTE` gathers RNA genes,
/// satellites and other non-TE references.
pub const FEATURE_CLASSES: &[&str] = &["DNA", "RC", "LTR", "LINE", "SINE", "PLE", "Retroposon"];

/// Domain roles reported as presence columns.
const ROLES: &[&str] = &[
    "GAG", "PR", "RT", "RH", "INT", "CHD", "APE", "GIY", "YR", "POLB", "HEL", "HELI", "HUH",
    "TPASE", "ULP",
];

/// Evidence for one family, as gathered by the classify pipeline.
pub struct FeatureInputs<'a> {
    pub family_id: &'a str,
    pub original_label: &'a str,
    pub label: &'a str,
    pub disposition: &'a str,
    pub confidence: &'a str,
    pub blastx: &'a HomologySummary,
    pub hmm: &'a HomologySummary,
    pub nucleotide: &'a HomologySummary,
    pub domains: &'a DomainEvidence,
    pub structure: &'a StructuralEvidence,
    pub copies: Option<&'a CopyEvidence>,
    pub host_domains: &'a [(String, f64)],
    /// A nucleotide hit lies in an element embedded in its reference (mosaics).
    pub embedded_hit: bool,
}

fn role_code(role: Role) -> &'static str {
    match role {
        Role::Gag => "GAG",
        Role::Protease | Role::PaoProtease => "PR",
        Role::ReverseTranscriptase => "RT",
        Role::RnaseH => "RH",
        Role::Integrase => "INT",
        Role::Chromo => "CHD",
        Role::Endonuclease => "APE",
        Role::GiyYig => "GIY",
        Role::TyrosineRecombinase => "YR",
        Role::PolB => "POLB",
        Role::HelitronHelicase => "HEL",
        Role::Helicase => "HELI",
        Role::Huh => "HUH",
        Role::Transposase(_) => "TPASE",
        Role::UlpProtease => "ULP",
    }
}

/// Feature group of a label: one of [`FEATURE_CLASSES`], `nonTE`, or none.
fn group(label: &str) -> Option<&'static str> {
    if is_non_te(label) {
        return Some("nonTE");
    }
    let class = canonical_class(label);
    FEATURE_CLASSES.iter().copied().find(|c| *c == class)
}

/// Best bits, consensus coverage and identity per class group.
fn by_class(summary: &HomologySummary, class: &str) -> (f64, f64, f64) {
    summary
        .classes
        .iter()
        .filter(|c| group(&c.class) == Some(class))
        .fold((0.0, 0.0, 0.0), |(b, cov, id), c| {
            (
                b.max(c.best_bitscore),
                cov.max(c.query_coverage),
                if c.best_bitscore > b {
                    c.best_identity
                } else {
                    id
                },
            )
        })
}

pub fn header() -> String {
    let mut h = vec![
        "family_id",
        "original_label",
        "label",
        "disposition",
        "confidence",
    ]
    .into_iter()
    .map(str::to_owned)
    .collect::<Vec<_>>();
    for class in FEATURE_CLASSES {
        h.extend([
            format!("prot_bits_{class}"),
            format!("prot_cov_{class}"),
            format!("hmm_bits_{class}"),
        ]);
    }
    for class in FEATURE_CLASSES.iter().chain(&["nonTE"]) {
        h.extend([
            format!("nucl_bits_{class}"),
            format!("nucl_cov_{class}"),
            format!("nucl_id_{class}"),
        ]);
    }
    h.extend(
        [
            "prot_unclassified_bits",
            "nucl_embedded_hit",
            "prot_chimeric",
            "nucl_chimeric",
        ]
        .map(str::to_owned),
    );
    for class in FEATURE_CLASSES {
        h.push(format!("dom_class_{class}"));
    }
    h.push("dom_novel".into());
    h.extend(ROLES.iter().map(|r| format!("dom_{r}")));
    h.extend(
        [
            "host_domains",
            "host_best_score",
            "length",
            "gc",
            "ltr_length",
            "ltr_identity",
            "tir_length",
            "tir_identity",
            "polya_3",
            "polyt_5",
            "tandem_period",
            "tandem_identity",
            "tandem_copies",
            "longest_orf_aa",
            "n_orfs",
            "tg_ca_termini",
            "helitron_termini",
            "has_copies",
            "copy_loci",
            "copy_full_length",
            "copy_identity",
            "genome_bp",
            "tsd_length",
            "tsd_support",
            "tsd_null",
            "flank_similarity",
            "helitron_ends",
            "helitron_ends_null",
        ]
        .map(str::to_owned),
    );
    h.join("\t")
}

pub fn row(f: &FeatureInputs) -> String {
    let mut cols: Vec<String> = vec![
        f.family_id.into(),
        f.original_label.into(),
        f.label.into(),
        f.disposition.into(),
        f.confidence.into(),
    ];
    let num = |v: f64| {
        if v == 0.0 {
            "0".to_owned()
        } else {
            format!("{v:.4}")
                .trim_end_matches('0')
                .trim_end_matches('.')
                .to_owned()
        }
    };
    let flag = |b: bool| if b { "1" } else { "0" }.to_owned();
    let opt = |v: Option<f64>| v.map(num).unwrap_or_default();
    for class in FEATURE_CLASSES {
        let (bits, cov, _) = by_class(f.blastx, class);
        let (hmm, _, _) = by_class(f.hmm, class);
        cols.extend([num(bits), num(cov), num(hmm)]);
    }
    for class in FEATURE_CLASSES.iter().chain(&["nonTE"]) {
        let (bits, cov, id) = by_class(f.nucleotide, class);
        cols.extend([num(bits), num(cov), num(id)]);
    }
    cols.extend([
        num(f.blastx.best_unclassified_bitscore),
        flag(f.embedded_hit),
        flag(!f.blastx.chimeric_classes.is_empty()),
        flag(!f.nucleotide.chimeric_classes.is_empty()),
    ]);
    let inferred = f.domains.inferred.as_deref().and_then(group);
    for class in FEATURE_CLASSES {
        cols.push(flag(inferred == Some(*class)));
    }
    cols.push(flag(f.domains.novel.is_some()));
    let roles = f
        .domains
        .hits
        .iter()
        .map(|h| role_code(h.role))
        .collect::<BTreeSet<_>>();
    cols.extend(ROLES.iter().map(|r| flag(roles.contains(r))));

    let s = f.structure;
    let c = f.copies;
    cols.extend([
        f.host_domains.len().to_string(),
        num(f.host_domains.iter().map(|(_, v)| *v).fold(0.0, f64::max)),
        s.length.to_string(),
        num(s.composition.gc),
        opt(s.ltr.as_ref().map(|r| r.length as f64)),
        opt(s.ltr.as_ref().map(|r| r.identity)),
        opt(s.tir.as_ref().map(|r| r.length as f64)),
        opt(s.tir.as_ref().map(|r| r.identity)),
        flag(s.tail == Some(Tail::PolyA3)),
        flag(s.tail == Some(Tail::PolyT5)),
        opt(s.tandem.as_ref().map(|t| t.period as f64)),
        opt(s.tandem.as_ref().map(|t| t.identity)),
        opt(s.tandem.as_ref().map(|t| t.copies)),
        s.longest_orf_aa().to_string(),
        s.orfs.len().to_string(),
        flag(s.tg_ca_termini),
        flag(s.helitron_termini),
        flag(c.is_some()),
        opt(c.map(|c| c.loci as f64)),
        opt(c.map(|c| c.full_length as f64)),
        opt(c.map(|c| c.mean_identity)),
        opt(c.map(|c| c.genome_bp as f64)),
        opt(c.and_then(|c| c.tsd.as_ref()).map(|t| t.length as f64)),
        opt(c.and_then(|c| c.tsd.as_ref()).map(|t| t.support)),
        opt(c.and_then(|c| c.tsd.as_ref()).map(|t| t.null_support)),
        opt(c.and_then(|c| c.flank_similarity)),
        opt(c.and_then(|c| c.helitron_ends.as_ref()).map(|h| h.support)),
        opt(c
            .and_then(|c| c.helitron_ends.as_ref())
            .map(|h| h.null_support)),
    ]);
    let mut line = String::new();
    write!(line, "{}", cols.join("\t")).expect("writing to a String cannot fail");
    line
}
