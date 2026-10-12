//! RepeatModeler / RepeatMasker library parsing.
//!
//! Library identifiers have the form `name#Class/Superfamily`, for example
//! `rnd-1_family-12#Unknown` or `ltr-1_family-3#LTR/Gypsy`.

use std::{
    collections::{HashMap, HashSet},
    path::Path,
};

use bio::io::fasta;

use crate::Result;

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

/// Read a library. Identifiers that occur more than once (seen in libraries
/// concatenated from several runs) are made unique by suffixing the name with
/// `_dup<n>`, keeping the label; identical duplicates are dropped.
pub fn read_library(path: &Path) -> Result<Vec<LibraryFamily>> {
    let reader = fasta::Reader::from_file(path)?;
    let mut seen: HashMap<String, Vec<Vec<u8>>> = HashMap::new();
    let mut ids = HashSet::new();
    let (mut renamed, mut dropped) = (0usize, 0usize);
    let mut families = Vec::new();
    for record in reader.records() {
        let record = record?;
        let original = record.id().to_owned();
        let seq = record.seq().to_ascii_uppercase();
        let copies = seen.entry(original.clone()).or_default();
        if copies.contains(&seq) {
            dropped += 1;
            continue;
        }
        copies.push(seq.clone());
        let mut id = original.clone();
        if copies.len() > 1 {
            let (name, label) = split_label(&original);
            let mut n = copies.len();
            loop {
                id = match &label {
                    Some(l) => format!("{name}_dup{n}#{l}"),
                    None => format!("{name}_dup{n}"),
                };
                if !ids.contains(&id) {
                    break;
                }
                n += 1;
            }
            renamed += 1;
        }
        ids.insert(id.clone());
        let (name, label) = split_label(&id);
        families.push(LibraryFamily {
            id,
            name,
            label,
            description: record.desc().map(str::to_owned),
            seq,
        });
    }
    if renamed + dropped > 0 {
        eprintln!(
            "Warning: {}: {renamed} repeated identifiers renamed with _dup<n>, {dropped} identical duplicates dropped",
            path.display()
        );
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
    fn renames_repeated_identifiers_and_drops_identical_duplicates() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            f,
            ">a#LINE/L2\nACGT\n>a#LINE/L2\nACGT\n>a#LINE/L2\nTTTT\n>a#LINE/L2\nGGGG\n>b\nAC"
        )
        .unwrap();
        let lib = read_library(f.path()).unwrap();
        let ids = lib.iter().map(|f| f.id.as_str()).collect::<Vec<_>>();
        assert_eq!(ids, ["a#LINE/L2", "a_dup2#LINE/L2", "a_dup3#LINE/L2", "b"]);
        assert_eq!(lib[1].label.as_deref(), Some("LINE/L2"));
        assert_eq!(lib[1].name, "a_dup2");
    }
}
