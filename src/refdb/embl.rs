//! Streaming parser for Repbase EMBL releases.
//!
//! Only the fields the reference database needs are kept: the identifier,
//! keywords (which carry the Repbase classification), species, nucleotide
//! sequence and any CDS protein translations.

use std::io::BufRead;

use crate::{Error, ErrorKind, Result};

#[derive(Clone, Debug, Default, PartialEq)]
pub struct EmblEntry {
    pub id: String,
    /// `KW` lines joined with spaces, e.g. `Gypsy; LTR Retrotransposon; ...`.
    pub keywords: String,
    pub species: Option<String>,
    pub seq: Vec<u8>,
    pub translations: Vec<String>,
}

/// Parse every entry from an EMBL stream.
///
/// Real releases contain small defects (stray bytes, doubled `//`
/// terminators), so invalid UTF-8 is replaced and empty records are skipped.
pub fn parse_embl(mut reader: impl BufRead) -> Result<Vec<EmblEntry>> {
    let mut entries = Vec::new();
    let mut entry = EmblEntry::default();
    let mut in_sequence = false;
    // an open `/translation="...` qualifier that continues over several lines
    let mut translation: Option<String> = None;
    let mut buffer = Vec::new();
    loop {
        buffer.clear();
        if reader.read_until(b'\n', &mut buffer)? == 0 {
            break;
        }
        let line = String::from_utf8_lossy(&buffer);
        let line = line.trim_end_matches(['\n', '\r']);
        if line.starts_with("//") {
            if !entry.id.is_empty() {
                entries.push(std::mem::take(&mut entry));
            }
            entry = EmblEntry::default();
            in_sequence = false;
            translation = None;
            continue;
        }
        if in_sequence {
            entry.seq.extend(
                line.bytes()
                    .filter(u8::is_ascii_alphabetic)
                    .map(|b| b.to_ascii_uppercase()),
            );
            continue;
        }
        let (tag, value) = (
            line.get(..2).unwrap_or(""),
            line.get(5..).unwrap_or("").trim(),
        );
        match tag {
            "ID" => entry.id = value.split_whitespace().next().unwrap_or("").to_owned(),
            "KW" => {
                if !entry.keywords.is_empty() {
                    entry.keywords.push(' ');
                }
                entry.keywords.push_str(value);
            }
            "OS" if entry.species.is_none() => entry.species = Some(value.to_owned()),
            "SQ" => in_sequence = true,
            "FT" => {
                // qualifier indentation varies between releases, so do not rely on columns
                let qualifier = line.get(2..).unwrap_or("").trim();
                if let Some(open) = translation.as_mut() {
                    open.push_str(qualifier.trim_end_matches('"'));
                    if qualifier.ends_with('"') {
                        entry.translations.extend(translation.take());
                    }
                } else if let Some(rest) = qualifier.strip_prefix("/translation=\"") {
                    let protein = rest.trim_end_matches('"').to_owned();
                    if rest.ends_with('"') {
                        entry.translations.push(protein);
                    } else {
                        translation = Some(protein);
                    }
                }
            }
            _ => {}
        }
    }
    if !entry.id.is_empty() {
        return Err(Error::new(ErrorKind::GenericCli(format!(
            "unterminated EMBL record {}",
            entry.id
        ))));
    }
    Ok(entries)
}

#[cfg(test)]
mod tests {
    use super::*;

    const RECORD: &str = "ID   Gypsy-1_XX   repbase; DNA; INV; 30 BP.
XX
KW   Gypsy; LTR Retrotransposon; Transposable Element;
KW   Gypsy-1_XX.
XX
OS   Testus exampli
FT   CDS             1..30
FT                   /product=\"Gypsy-1_XX_1p\"
FT                 /translation=\"MKVLA
FT                   WWHHK\"
FT   CDS             1..9
FT                   /translation=\"MKV\"
XX
SQ   Sequence 30 BP; 8 A; 7 C; 8 G; 7 T; 0 other;
     acgtacgtac gtacgtacgt               20
     acgtacgtac                          30
//
";

    #[test]
    fn parses_ids_keywords_species_sequence_and_translations() {
        let entries = parse_embl(RECORD.repeat(2).as_bytes()).unwrap();
        assert_eq!(entries.len(), 2);
        let e = &entries[0];
        assert_eq!(e.id, "Gypsy-1_XX");
        assert_eq!(
            e.keywords,
            "Gypsy; LTR Retrotransposon; Transposable Element; Gypsy-1_XX."
        );
        assert_eq!(e.species.as_deref(), Some("Testus exampli"));
        assert_eq!(e.seq.len(), 30);
        assert!(e.seq.starts_with(b"ACGTACGT"));
        assert_eq!(e.translations, vec!["MKVLAWWHHK", "MKV"]);
    }

    #[test]
    fn tolerates_doubled_terminators_and_invalid_utf8() {
        let mut data = RECORD.as_bytes().to_vec();
        data.splice(data.len() - 3..data.len() - 3, [0xffu8, b'\n']);
        data.extend_from_slice(b"//\n");
        data.extend_from_slice(RECORD.as_bytes());
        let entries = parse_embl(data.as_slice()).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].seq.len(), 30);
    }

    #[test]
    fn rejects_unterminated_records() {
        assert!(parse_embl("ID   X   repbase;\nSQ\n acgt\n".as_bytes()).is_err());
    }
}
