//! RepeatModeler / RepeatMasker library parsing.
//!
//! Library identifiers have the form `name#Class/Superfamily`, for example
//! `rnd-1_family-12#Unknown` or `ltr-1_family-3#LTR/Gypsy`.

use std::{collections::HashSet, path::Path};

use bio::io::fasta;

use crate::{Error, ErrorKind, Result};

#[derive(Clone, Debug, PartialEq)]
pub struct LibraryFamily {
    /// Full FASTA identifier, including any `#label`.
    pub id: String,
    /// Identifier without the label.
    pub name: String,
    /// RepeatMasker label after `#`, if present.
    pub label: Option<String>,
    pub description: Option<String>,
    pub seq: Vec<u8>,
}

impl LibraryFamily {
    /// True when the family has no label or a label that starts with `Unknown`
    /// (`Unknown`, `Unknown/Unknown`, ...).
    pub fn is_unknown(&self) -> bool {
        self.label
            .as_deref()
            .is_none_or(|l| l.eq_ignore_ascii_case("unknown") || l.starts_with("Unknown"))
    }
}

/// Split `name#Class/Superfamily` into the name and label.
pub fn split_label(id: &str) -> (String, Option<String>) {
    match id.split_once('#') {
        Some((name, label)) if !label.is_empty() => (name.to_owned(), Some(label.to_owned())),
        Some((name, _)) => (name.to_owned(), None),
        None => (id.to_owned(), None),
    }
}

pub fn read_library(path: &Path) -> Result<Vec<LibraryFamily>> {
    let reader = fasta::Reader::from_file(path)?;
    let mut seen = HashSet::new();
    let mut families = Vec::new();
    for record in reader.records() {
        let record = record?;
        let id = record.id().to_owned();
        if !seen.insert(id.clone()) {
            return Err(Error::new(ErrorKind::GenericCli(format!(
                "duplicate library identifier {id} in {}",
                path.display()
            ))));
        }
        let (name, label) = split_label(&id);
        families.push(LibraryFamily {
            id,
            name,
            label,
            description: record.desc().map(str::to_owned),
            seq: record.seq().to_ascii_uppercase(),
        });
    }
    Ok(families)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn splits_repeatmodeler_labels() {
        assert_eq!(
            split_label("rnd-1_family-12#Unknown"),
            ("rnd-1_family-12".into(), Some("Unknown".into()))
        );
        assert_eq!(
            split_label("ltr-1_family-3#LTR/Gypsy"),
            ("ltr-1_family-3".into(), Some("LTR/Gypsy".into()))
        );
        assert_eq!(split_label("plain"), ("plain".into(), None));
    }

    #[test]
    fn reads_library_and_flags_unknowns() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            f,
            ">rnd-1_family-0#Unknown ( RepeatScout Family Size = 40 )\nacgt\n>rnd-1_family-1#LINE/L1\nACGT\n>nolabel\nAC"
        )
        .unwrap();
        let lib = read_library(f.path()).unwrap();
        assert_eq!(lib.len(), 3);
        assert!(lib[0].is_unknown());
        assert_eq!(lib[0].seq, b"ACGT");
        assert!(lib[0]
            .description
            .as_deref()
            .unwrap()
            .contains("RepeatScout"));
        assert!(!lib[1].is_unknown());
        assert!(lib[2].is_unknown());
    }

    #[test]
    fn rejects_duplicate_identifiers() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, ">a#Unknown\nACGT\n>a#Unknown\nACGT").unwrap();
        assert!(read_library(f.path()).is_err());
    }
}
