//! Mapping TE protein names to RepeatMasker-style `Class/Superfamily` labels.
//!
//! Most RepBase protein records carry only a name (`Gypsy-38_AnFu-I_1p`), so the
//! class is inferred from an ordered rule table. Records that already carry a
//! RepeatMasker label (`RTE-1_Ec32_EN-RT#LINE/RTE`) use that label directly.
//! Names that no rule matches return `None`; callers must report them rather
//! than guess.

/// How a rule keyword is matched against the upper-cased protein name.
#[derive(Clone, Copy)]
enum Match {
    /// The name starts with the keyword, followed by a non-letter or the end.
    Prefix,
    /// The keyword occurs anywhere in the name.
    Contains,
}

use Match::{Contains, Prefix};

/// Ordered rules: the first match wins, so specific rules precede generic ones.
const RULES: &[(&str, Match, &str)] = &[
    // --- LTR retrotransposons ---
    ("ATHILA", Contains, "LTR/Gypsy"),
    ("ATGP", Prefix, "LTR/Gypsy"),
    ("ATLANTYS", Prefix, "LTR/Gypsy"),
    ("GYPSY", Contains, "LTR/Gypsy"),
    ("GYPS", Prefix, "LTR/Gypsy"),
    ("GYP", Prefix, "LTR/Gypsy"),
    ("TY3", Prefix, "LTR/Gypsy"),
    ("OGRE", Contains, "LTR/Gypsy"),
    ("RIRE", Prefix, "LTR/Gypsy"),
    ("COPIA", Contains, "LTR/Copia"),
    ("COP", Prefix, "LTR/Copia"),
    ("SHACOP", Prefix, "LTR/Copia"),
    ("TY1", Prefix, "LTR/Copia"),
    ("RLG", Prefix, "LTR/Gypsy"),
    ("RLC", Prefix, "LTR/Copia"),
    ("BEL", Prefix, "LTR/Pao"),
    ("PAO", Prefix, "LTR/Pao"),
    ("ERV", Contains, "LTR/ERV"),
    ("IAP", Prefix, "LTR/ERVK"),
    ("MMTV", Prefix, "LTR/ERVK"),
    ("SFV", Prefix, "LTR/ERV"),
    ("SRV", Prefix, "LTR/ERV"),
    ("RSV", Prefix, "LTR/ERV"),
    ("CAULIMO", Contains, "LTR/Caulimovirus"),
    ("CAULIV", Contains, "LTR/Caulimovirus"),
    ("DIRS", Contains, "LTR/DIRS"),
    ("NGARO", Contains, "LTR/Ngaro"),
    ("VIPER", Contains, "LTR/DIRS"),
    ("PAT", Prefix, "LTR/DIRS"),
    // --- Penelope-like ---
    ("PENELOPE", Contains, "PLE/Penelope"),
    ("POSEIDON", Prefix, "PLE/Penelope"),
    ("NEPTUNE", Prefix, "PLE/Penelope"),
    ("ATHENA", Prefix, "PLE/Penelope"),
    ("NAIAD", Prefix, "PLE/Penelope"),
    ("CHLAMYS", Prefix, "PLE/Penelope"),
    // --- non-LTR (LINE) retrotransposons ---
    ("ATLINE", Prefix, "LINE"),
    ("LINE1", Prefix, "LINE/L1"),
    ("L1", Prefix, "LINE/L1"),
    ("TX1", Prefix, "LINE/L1"),
    ("L2", Prefix, "LINE/L2"),
    ("CR1", Prefix, "LINE/CR1"),
    ("REX", Prefix, "LINE/Rex-Babar"),
    ("RTEX", Prefix, "LINE/RTE-X"),
    ("RTE", Prefix, "LINE/RTE"),
    ("BOVB", Prefix, "LINE/RTE-BovB"),
    ("JOCKEY", Prefix, "LINE/I-Jockey"),
    ("JOCK", Prefix, "LINE/I-Jockey"),
    ("I", Prefix, "LINE/I"),
    ("R1", Prefix, "LINE/R1"),
    ("R2", Prefix, "LINE/R2"),
    ("R4", Prefix, "LINE/R4"),
    ("NESL", Prefix, "LINE/R2-NeSL"),
    ("CRE", Prefix, "LINE/CRE"),
    ("TAD", Prefix, "LINE/Tad1"),
    ("LOA", Prefix, "LINE/I"),
    ("OUTCAST", Prefix, "LINE/Tad1"),
    ("VINGI", Prefix, "LINE/R2-Vingi"),
    ("INGI", Prefix, "LINE/I"),
    ("DONG", Prefix, "LINE/R4-Dong"),
    ("PROTO", Prefix, "LINE/Proto"),
    ("KIRI", Prefix, "LINE/Kiri"),
    ("NIMB", Prefix, "LINE/Nimb"),
    ("NIMBX", Prefix, "LINE/Nimb"),
    ("DAPHNE", Prefix, "LINE/Daphne"),
    ("CRACK", Prefix, "LINE/Crack"),
    ("HERO", Prefix, "LINE/Hero"),
    ("TART", Prefix, "LINE/I-Jockey"),
    ("SART", Prefix, "LINE/R1"),
    ("TRAS", Prefix, "LINE/R1"),
    ("LINE", Prefix, "LINE"),
    ("NONLTR", Prefix, "LINE"),
    // --- rolling-circle and self-synthesising DNA elements ---
    ("HELITRON", Contains, "RC/Helitron"),
    ("ATREP", Prefix, "RC/Helitron"),
    ("POLINTON", Contains, "DNA/Maverick"),
    ("MAVERICK", Contains, "DNA/Maverick"),
    ("CRYPTON", Contains, "DNA/Crypton"),
    // --- DDE DNA transposons ---
    ("PIGGYB", Contains, "DNA/PiggyBac"),
    ("HARBINGER", Contains, "DNA/PIF-Harbinger"),
    ("HARB", Prefix, "DNA/PIF-Harbinger"),
    ("MARINER", Contains, "DNA/TcMar-Mariner"),
    ("TCMAR", Prefix, "DNA/TcMar"),
    ("TC1", Prefix, "DNA/TcMar-Tc1"),
    ("TC2", Prefix, "DNA/TcMar-Tc2"),
    ("POGO", Prefix, "DNA/TcMar-Pogo"),
    ("MINOS", Prefix, "DNA/TcMar-Tc1"),
    ("ENSPM", Contains, "DNA/CMC-EnSpm"),
    ("CACTA", Contains, "DNA/CMC-EnSpm"),
    ("CHAPAEV", Contains, "DNA/CMC-Chapaev"),
    ("TRANSIB", Contains, "DNA/CMC-Transib"),
    ("MUDR", Contains, "DNA/MULE-MuDR"),
    ("MUTATOR", Contains, "DNA/MULE"),
    ("ATMU", Prefix, "DNA/MULE"),
    ("ATHAT", Prefix, "DNA/hAT"),
    ("HATM", Prefix, "DNA/hAT"),
    ("HATX", Prefix, "DNA/hAT"),
    ("HATW", Prefix, "DNA/hAT"),
    ("HATV", Prefix, "DNA/hAT"),
    ("ZHAT", Prefix, "DNA/hAT"),
    ("VIHAT", Prefix, "DNA/hAT"),
    ("HAT", Prefix, "DNA/hAT"),
    ("SPIN", Prefix, "DNA/hAT"),
    ("ACADEM", Contains, "DNA/Academ"),
    ("KOLOBOK", Contains, "DNA/Kolobok"),
    ("GINGER", Prefix, "DNA/Ginger"),
    ("SOLA", Prefix, "DNA/Sola"),
    ("ZATOR", Prefix, "DNA/Zator"),
    ("ZISUPTON", Prefix, "DNA/Zisupton"),
    ("NOVOSIB", Prefix, "DNA/Novosib"),
    ("MERLIN", Prefix, "DNA/Merlin"),
    ("DADA", Prefix, "DNA/Dada"),
    ("PLAVAKA", Prefix, "DNA/Plavaka"),
    ("ISL", Prefix, "DNA/IS3EU"),
    ("DNA", Prefix, "DNA"),
];

/// Infer a RepeatMasker-style class label from a TE protein record name.
///
/// A `#Class/Family` suffix takes precedence; otherwise the ordered rule table
/// is applied to the upper-cased name.
pub fn class_from_protein_name(name: &str) -> Option<String> {
    if let Some((_, label)) = name.split_once('#') {
        let label = label.split_whitespace().next().unwrap_or("");
        if !label.is_empty() && !label.starts_with("Unknown") {
            return Some(label.to_owned());
        }
    }
    // Some RepBase names carry a species-group prefix (`Ag-CR1-22_1p`,
    // `Adg_Copia-70_AA-I_1p`), so the rules are also tried after it.
    let upper = name.to_ascii_uppercase();
    let mut candidates = vec![upper.as_str()];
    if let Some(index) = upper.find(['-', '_']) {
        if index <= 4 {
            candidates.push(&upper[index + 1..]);
        }
    }
    for candidate in candidates {
        for &(keyword, mode, class) in RULES {
            let hit = match mode {
                Contains => candidate.contains(keyword),
                Prefix => {
                    candidate.starts_with(keyword)
                        && !candidate[keyword.len()..]
                            .chars()
                            .next()
                            .is_some_and(|c| c.is_ascii_alphabetic())
                }
            };
            if hit {
                return Some(class.to_owned());
            }
        }
    }
    None
}

/// The top-level class of a label, e.g. `LTR` for `LTR/Gypsy`.
pub fn top_level(label: &str) -> &str {
    label.split('/').next().unwrap_or(label)
}

/// [`top_level`] with legacy class names mapped to current ones:
/// RepeatMasker formerly filed Penelope-like elements under `LINE`.
pub fn canonical_class(label: &str) -> &str {
    if label.starts_with("LINE/Penelope") {
        "PLE"
    } else {
        top_level(label)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn explicit_labels_take_precedence() {
        assert_eq!(
            class_from_protein_name("RTEX-1_Ec32_EN-RT#LINE/RTE-X").as_deref(),
            Some("LINE/RTE-X")
        );
    }

    #[test]
    fn resolves_common_repbase_names() {
        for (name, class) in [
            ("Gypsy-38_AnFu-I_1p", "LTR/Gypsy"),
            ("GYPSY27_AGp", "LTR/Gypsy"),
            ("COPIA5_AGp", "LTR/Copia"),
            ("BEL-12_AnFu-I_1p", "LTR/Pao"),
            ("CR1-6_AG-ORF2p", "LINE/CR1"),
            ("Ag-CR1-22_1p", "LINE/CR1"),
            ("Adg_Copia-70_AA-I_1p", "LTR/Copia"),
            ("LINE1-6_OS_1p", "LINE/L1"),
            ("L1-5_xx", "LINE/L1"),
            ("I-1_DR", "LINE/I"),
            ("RTEX-2_xx", "LINE/RTE-X"),
            ("hATm-1_xx", "DNA/hAT"),
            ("hAT-5_xx", "DNA/hAT"),
            ("Mariner-3_xx", "DNA/TcMar-Mariner"),
            ("piggyBacX-1_xx", "DNA/PiggyBac"),
            ("Helitron-1_xx", "RC/Helitron"),
            ("Polinton-2_xx", "DNA/Maverick"),
            ("MicrocebusERV1", "LTR/ERV"),
            ("DNA-TTAA-3_NV_1p", "DNA"),
        ] {
            assert_eq!(
                class_from_protein_name(name).as_deref(),
                Some(class),
                "{name}"
            );
        }
    }

    #[test]
    fn prefix_rules_respect_word_boundaries() {
        // `I` must not match `Itch`, and `L1` must not match `L1x`-like words
        assert_eq!(class_from_protein_name("Itch-1_xx"), None);
        assert_eq!(class_from_protein_name("Utopia-3_xx"), None);
    }

    #[test]
    fn legacy_penelope_is_ple() {
        assert_eq!(canonical_class("LINE/Penelope"), "PLE");
        assert_eq!(canonical_class("LINE/L2"), "LINE");
        assert_eq!(
            class_from_protein_name("Penelope-1_DV_1p").as_deref(),
            Some("PLE/Penelope")
        );
    }

    #[test]
    fn top_level_class() {
        assert_eq!(top_level("LTR/Gypsy"), "LTR");
        assert_eq!(top_level("DNA"), "DNA");
    }
}
