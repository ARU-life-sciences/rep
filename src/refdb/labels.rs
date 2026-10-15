//! Repbase EMBL keywords → RepeatMasker-style `Class/Superfamily` labels.
//!
//! Repbase `KW` lines start with `superfamily; class; ...`. Superfamilies with
//! an established RepeatMasker name map to it. Other superfamilies keep their
//! Repbase name under the mapped class (`LINE/Daphne`), so no information is
//! lost. Records that give neither a usable class nor superfamily are `None`.

use crate::classify::class_from_protein_name;

/// RepeatMasker names for Repbase superfamilies (matched case-insensitively).
const SUPERFAMILIES: &[(&str, &str)] = &[
    ("Gypsy", "LTR/Gypsy"),
    ("Copia", "LTR/Copia"),
    ("BEL", "LTR/Pao"),
    ("DIRS", "LTR/DIRS"),
    ("Ngaro", "LTR/Ngaro"),
    ("ERV1", "LTR/ERV1"),
    ("ERV2", "LTR/ERVK"),
    ("ERV3", "LTR/ERVL"),
    ("ERV4", "LTR/ERV4"),
    ("Caulimoviridae", "LTR/Caulimovirus"),
    ("hAT", "DNA/hAT"),
    ("Mariner/Tc1", "DNA/TcMar"),
    ("MuDR", "DNA/MULE-MuDR"),
    ("Harbinger", "DNA/PIF-Harbinger"),
    ("Helitron", "RC/Helitron"),
    ("EnSpm/CACTA", "DNA/CMC-EnSpm"),
    ("EnSpm", "DNA/CMC-EnSpm"),
    ("Chapaev", "DNA/CMC-Chapaev"),
    ("Transib", "DNA/CMC-Transib"),
    ("piggyBac", "DNA/PiggyBac"),
    ("Polinton", "DNA/Maverick"),
    ("Sola1", "DNA/Sola-1"),
    ("Sola2", "DNA/Sola-2"),
    ("Sola3", "DNA/Sola-3"),
    ("Ginger1", "DNA/Ginger-1"),
    ("Ginger2/TDD", "DNA/Ginger-2"),
    ("L1", "LINE/L1"),
    ("Tx1", "LINE/L1-Tx1"),
    ("L2", "LINE/L2"),
    ("L2B", "LINE/L2"),
    ("CR1", "LINE/CR1"),
    ("RTE", "LINE/RTE"),
    ("RTEX", "LINE/RTE-X"),
    ("Jockey", "LINE/I-Jockey"),
    ("I", "LINE/I"),
    ("Loa", "LINE/I"),
    ("Ingi", "LINE/I"),
    ("R2", "LINE/R2"),
    ("NeSL", "LINE/R2-NeSL"),
    ("Rex1", "LINE/Rex-Babar"),
    ("Penelope", "PLE/Penelope"),
    ("Penelope/Poseidon", "PLE/Penelope"),
    ("SINE2/tRNA", "SINE/tRNA"),
    ("SINE1/7SL", "SINE/7SL"),
    ("SINE3/5S", "SINE/5S"),
    ("SINEU/snRNA", "SINE/U"),
    ("SINE", "SINE"),
    ("SAT", "Satellite"),
    ("MSAT", "Satellite"),
    ("tRNA", "tRNA"),
    ("rRNA", "rRNA"),
    ("snRNA", "snRNA"),
];

/// Repbase class field → RepeatMasker class.
const CLASSES: &[(&str, &str)] = &[
    ("LTR Retrotransposon", "LTR"),
    ("Endogenous Retrovirus", "LTR"),
    ("Integrated Virus", "LTR"),
    ("Non-LTR Retrotransposon", "LINE"),
    ("SINE", "SINE"),
    ("DNA transposon", "DNA"),
    ("Sola", "DNA"),
    ("Crypton", "DNA"),
    ("Penelope", "PLE"),
    ("Satellite", "Satellite"),
    ("Simple Repeat", "Simple_repeat"),
    ("Multicopy gene", "Other"),
];

fn lookup<'a>(table: &'a [(&str, &str)], key: &str) -> Option<&'a str> {
    table
        .iter()
        .find(|(k, _)| k.eq_ignore_ascii_case(key))
        .map(|&(_, v)| v)
}

/// Map a Repbase `KW` string to a RepeatMasker label.
pub fn repbase_label(keywords: &str) -> Option<String> {
    let mut fields = keywords
        .split(';')
        .map(|f| f.trim().trim_end_matches('.').trim());
    let superfamily = fields.next().unwrap_or("");
    let class_field = fields.next().unwrap_or("");

    // `DNA transposon; Transposable Element`: the class sits in the first field
    if class_field.eq_ignore_ascii_case("Transposable Element") {
        return lookup(CLASSES, superfamily).map(str::to_owned);
    }
    let class = lookup(CLASSES, class_field);
    if class == Some("Simple_repeat") {
        return Some("Simple_repeat".into());
    }
    // Crypton and Sola are classes in Repbase but DNA superfamilies in RepeatMasker
    if let Some(prefix) = ["Crypton", "Sola"]
        .iter()
        .find(|p| class_field.eq_ignore_ascii_case(p))
    {
        if let Some(label) = lookup(SUPERFAMILIES, superfamily) {
            return Some(label.to_owned());
        }
        let suffix = superfamily.strip_prefix(prefix).unwrap_or("");
        return Some(if suffix.is_empty() {
            format!("DNA/{prefix}")
        } else {
            format!("DNA/{prefix}-{suffix}")
        });
    }
    if class_field.eq_ignore_ascii_case("Non-LTR Retrotransposon")
        && superfamily.eq_ignore_ascii_case("SINE")
    {
        return Some("SINE".into());
    }
    if let Some(label) = lookup(SUPERFAMILIES, superfamily) {
        // a known superfamily must agree with the stated class, if any
        let consistent = class.is_none_or(|c| {
            label.starts_with(c)
                || (c == "Other" && matches!(label, "tRNA" | "rRNA" | "snRNA"))
                || (c == "LINE" && label.starts_with("PLE"))
                // Repbase files Helitrons under the DNA transposon class
                || (c == "DNA" && label.starts_with("RC/"))
        });
        if consistent {
            return Some(label.to_owned());
        }
    }
    match class {
        Some("Other") => Some(format!("Other/{}", sanitize(superfamily))),
        Some("Satellite") => Some("Satellite".into()),
        Some(c) if !superfamily.is_empty() => Some(format!("{c}/{}", sanitize(superfamily))),
        Some(c) => Some(c.to_owned()),
        // no recognised class: fall back to the name-based TE rules
        None => class_from_protein_name(superfamily),
    }
}

/// Lower-case words of a Repbase description or keyword line (hyphens kept,
/// so `non-LTR` stays one word).
fn words(text: &str) -> Vec<String> {
    text.split(|c: char| !(c.is_ascii_alphanumeric() || c == '-'))
        .filter(|w| !w.is_empty())
        .map(str::to_ascii_lowercase)
        .collect()
}

/// Does the text name a transposable element or its proteins?
fn mentions_te(words: &[String]) -> bool {
    words.iter().any(|w| {
        w.contains("transpos")
            || w.contains("retro")
            || w.contains("viru")
            || w.contains("viral")
            || matches!(
                w.as_str(),
                "ltr"
                    | "line"
                    | "sine"
                    | "mite"
                    | "helitron"
                    | "element"
                    | "integrase"
                    | "transcriptase"
                    | "gag"
                    | "polyprotein"
                    | "non-ltr"
            )
    })
}

/// Why a Repbase entry is left out of the reference database, if it is.
///
/// Repbase also holds protein-coding multicopy genes (ribosomal proteins,
/// histones) and, in older entries, repeats described only by the host
/// gene or mRNA they were found in (`UBX_DM`: "Drosophila repeat, from
/// Ultrabithorax locus"; `PEN1`: "pen repeat mRNA"). As references they
/// make host genes look like repeats. Entries with a usable label are kept,
/// as are RNA genes and anything whose description names a TE.
pub fn repbase_exclusion(
    keywords: &str,
    description: &str,
    label: Option<&str>,
) -> Option<&'static str> {
    let kw = words(keywords);
    let de = words(description);
    let rna = |w: &[String]| {
        w.iter()
            .any(|x| x.ends_with("rna") || matches!(x.as_str(), "5s" | "7sl" | "srp" | "ncrna"))
    };
    let multicopy = keywords.split(';').any(|f| {
        f.trim()
            .trim_end_matches('.')
            .eq_ignore_ascii_case("Multicopy gene")
    });
    if multicopy && !matches!(label, Some("tRNA" | "rRNA" | "snRNA")) && !rna(&kw) && !rna(&de) {
        return Some("multicopy_gene");
    }
    if label.is_some() {
        return None;
    }
    let all = kw.iter().chain(&de).cloned().collect::<Vec<_>>();
    let host = de.iter().any(|w| {
        matches!(
            w.as_str(),
            "gene"
                | "genes"
                | "mrna"
                | "cdna"
                | "cds"
                | "locus"
                | "allele"
                | "mutant"
                | "intron"
                | "exon"
                | "promoter"
                | "utr"
                | "flanking"
                | "non-coding"
        )
    });
    (host && !mentions_te(&all)).then_some("host_gene_context")
}

/// A label from the free text of entries whose keywords [`repbase_label`]
/// cannot read: the old-style `LTR; retrotransposon; ...` keywords of the
/// appendix, `Simple repeat.`, RNA genes, and descriptions such as "Drosophila
/// retrotransposon 412 genome" with `LTR` among the keywords. Ambiguous text
/// (a bare "retrotransposon") stays unlabelled.
pub fn repbase_legacy_label(keywords: &str, description: &str) -> Option<String> {
    let w = words(&format!("{keywords} {description}"));
    let has = |x: &str| w.iter().any(|y| y == x);
    let text = w.join(" ");
    let label = if text.contains("simple repeat") {
        "Simple_repeat"
    } else if has("satellite") || has("minisatellite") {
        "Satellite"
    } else if has("helitron") {
        "RC/Helitron"
    } else if has("sine") {
        "SINE"
    } else if has("non-ltr") || has("line") || has("line-1") {
        "LINE"
    } else if has("ltr") || text.contains("endogenous retrovirus") {
        "LTR"
    } else if has("mite")
        || has("transposase")
        || text.contains("dna transposon")
        || (w.iter().any(|x| x == "transposon" || x == "transposons")
            && !w.iter().any(|x| x.contains("retro")))
    {
        "DNA"
    } else if has("trna") {
        "tRNA"
    } else if ["rrna", "5s", "18s", "28s"].iter().any(|x| has(x)) {
        "rRNA"
    } else if w
        .iter()
        .any(|x| x == "snrna" || (x.starts_with('u') && x.ends_with("-rna")))
    {
        "snRNA"
    } else if w
        .iter()
        .any(|x| x.ends_with("rna") && !matches!(x.as_str(), "mrna" | "cdna"))
    {
        "RNA"
    } else {
        return None;
    };
    Some(label.to_owned())
}

/// Make a Repbase superfamily name usable inside a RepeatMasker label.
fn sanitize(name: &str) -> String {
    name.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' {
                c
            } else {
                '-'
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn maps_common_repbase_keywords() {
        for (kw, label) in [
            (
                "Gypsy; LTR Retrotransposon; Transposable Element; X.",
                "LTR/Gypsy",
            ),
            ("BEL; LTR Retrotransposon;", "LTR/Pao"),
            ("Mariner/Tc1; DNA transposon;", "DNA/TcMar"),
            ("ERV2; Endogenous Retrovirus;", "LTR/ERVK"),
            ("DNA transposon; Transposable Element;", "DNA"),
            ("Non-LTR Retrotransposon; Transposable Element;", "LINE"),
            ("SINE2/tRNA; SINE;", "SINE/tRNA"),
            ("SINE; Non-LTR Retrotransposon;", "SINE"),
            ("Daphne; Non-LTR Retrotransposon;", "LINE/Daphne"),
            ("Penelope; Non-LTR Retrotransposon;", "PLE/Penelope"),
            ("Neptune; Penelope;", "PLE/Neptune"),
            ("Naiad/Chlamys; Penelope;", "PLE/Naiad-Chlamys"),
            ("Sola2; Sola;", "DNA/Sola-2"),
            ("CryptonV; Crypton;", "DNA/Crypton-V"),
            ("SAT; Satellite;", "Satellite"),
            ("Satellite; Simple Repeat;", "Simple_repeat"),
            ("tRNA; Multicopy gene;", "tRNA"),
            ("Caulimoviridae; Integrated Virus;", "LTR/Caulimovirus"),
            ("Troyka; LTR Retrotransposon;", "LTR/Troyka"),
            (
                "Helitron; DNA transposon; Transposable Element;",
                "RC/Helitron",
            ),
        ] {
            assert_eq!(repbase_label(kw).as_deref(), Some(label), "{kw}");
        }
    }

    #[test]
    fn legacy_entries_are_labelled_from_their_text() {
        for (kw, de, label) in [
            (
                "LTR; retrotransposon; reverse transcriptase; DM412.",
                "Drosophila retrotransposon 412 genome.",
                Some("LTR"),
            ),
            ("Simple repeat.", "(CA)n repeat", Some("Simple_repeat")),
            ("non-LTR retrotransposon; X.", "", Some("LINE")),
            ("MITE; X.", "", Some("DNA")),
            ("transposon; TZ86.", "maize transposon Tz86", Some("DNA")),
            ("retrotransposon; X.", "a retrotransposon", None),
            ("Repetitive sequence; UBX_DM.", "Drosophila repeat", None),
            (
                "Multicopy gene; ncRNA; RNA; RNA-U1_OS.",
                "U1-RNA gene sequence",
                Some("snRNA"),
            ),
            (
                "PEN1; Repetitive sequence.",
                "D.melanogaster pen repeat mRNA.",
                None,
            ),
        ] {
            assert_eq!(repbase_legacy_label(kw, de).as_deref(), label, "{kw} {de}");
        }
    }

    #[test]
    fn host_genes_are_excluded_but_tes_and_rna_genes_kept() {
        let ex = |kw, de, label| repbase_exclusion(kw, de, label);
        assert_eq!(
            ex(
                "Repetitive sequence; DMUBX; UBX_DM.",
                "Drosophila repeat, from Ultrabithorax locus (Ubx) of the Bithorax complex (BX-C).",
                None
            ),
            Some("host_gene_context")
        );
        assert_eq!(
            ex(
                "PEN1; Repetitive sequence.",
                "D.melanogaster pen repeat mRNA.",
                None
            ),
            Some("host_gene_context")
        );
        assert_eq!(
            ex(
                "Multicopy gene; 60S-L21_MM.",
                "Multiple copied gene: 60S ribosomal protein L21 gene.",
                None
            ),
            Some("multicopy_gene")
        );
        assert_eq!(
            ex(
                "Multicopy gene; ncRNA; RNA; RNA-U1_OS.",
                "U1-RNA gene sequence",
                Some("snRNA")
            ),
            None
        );
        assert_eq!(ex("tRNA; Multicopy gene;", "", Some("tRNA")), None);
        assert_eq!(
            ex(
                "R1B_DS.",
                "Dugesiella sp. retrotransposon R1 reverse transcriptase gene, partial cds.",
                None
            ),
            None
        );
        assert_eq!(
            ex("Gypsy; LTR Retrotransposon;", "a gene", Some("LTR/Gypsy")),
            None
        );
    }

    #[test]
    fn unusable_keywords_are_unmapped() {
        assert_eq!(repbase_label(""), None);
        assert_eq!(repbase_label("Nonautonomous; Transposable Element;"), None);
    }
}
