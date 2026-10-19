//! `rep relatives`: the families that `rep classify` runs called as TEs,
//! pooled per taxon, for `rep classify --relatives` on related genomes.
//!
//! A family is a reference when the cascade called it a TE (`known_te` or
//! `chimeric` with a TE label) or the learned classifier decided it; families
//! labelled from relatives or genomic neighbours are left out, so labels do not
//! travel along chains.

use std::{
    collections::{BTreeMap, HashMap},
    fs,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use bio::io::fasta;
use serde::Serialize;

use super::{classify_dirs, read_classification, read_taxa};
use crate::{classify::is_non_te, Result};

#[derive(Clone, Debug)]
pub struct RelativesConfig {
    /// Root directory searched for `rep classify` outputs.
    pub results: PathBuf,
    /// TSV with a species (`scientific_name` or `species`) column and `rank`.
    pub taxa: Option<PathBuf>,
    /// Taxonomy column that groups samples, e.g. `order`.
    pub rank: String,
    pub output_dir: PathBuf,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct RelativesSummary {
    pub samples: usize,
    /// Reference families written per taxon file.
    pub per_taxon: BTreeMap<String, usize>,
}

/// File stem for samples whose species has no taxon.
pub const UNASSIGNED: &str = "unassigned";

/// Is this `classification.tsv` row a reference for related genomes?
fn is_reference(row: &std::collections::HashMap<String, String>) -> bool {
    let get = |k: &str| row.get(k).map_or("", String::as_str);
    let label = get("label");
    match get("decided_by") {
        "model" => true,
        "relative" | "neighbour" => false,
        _ => {
            matches!(get("disposition"), "known_te" | "chimeric")
                && !label.is_empty()
                && label != "Unknown"
                && !is_non_te(label)
        }
    }
}

pub fn write_relatives(c: &RelativesConfig) -> Result<RelativesSummary> {
    fs::create_dir_all(&c.output_dir)?;
    let taxa = match &c.taxa {
        Some(p) => read_taxa(p, &c.rank)?,
        None => Default::default(),
    };
    let dirs = classify_dirs(&c.results);
    let mut files: BTreeMap<String, BufWriter<fs::File>> = BTreeMap::new();
    let mut summary = RelativesSummary {
        samples: dirs.len(),
        ..Default::default()
    };
    for dir in &dirs {
        let sample = dir
            .strip_prefix(&c.results)
            .unwrap_or(dir)
            .to_string_lossy()
            .into_owned();
        let species = sample.split('/').next().unwrap_or(&sample);
        let taxon = taxon_of(&taxa, species)
            .unwrap_or(UNASSIGNED)
            .replace(['/', ' '], "_");
        let rows = read_classification(&dir.join("classification.tsv"))?;
        let mut written = 0;
        for record in fasta::Reader::from_file(dir.join("queries.fa"))?.records() {
            let record = record?;
            let Some(row) = rows.get(record.id()).filter(|r| is_reference(r)) else {
                continue;
            };
            let name = record.id().split('#').next().unwrap_or(record.id());
            let out = match files.get_mut(&taxon) {
                Some(f) => f,
                None => files
                    .entry(taxon.clone())
                    .or_insert(BufWriter::new(fs::File::create(
                        c.output_dir.join(format!("{taxon}.fa")),
                    )?)),
            };
            writeln!(out, ">{}|{name}#{}", sample.replace('/', ":"), row["label"])?;
            out.write_all(&record.seq().to_ascii_uppercase())?;
            out.write_all(b"\n")?;
            written += 1;
        }
        *summary.per_taxon.entry(taxon).or_default() += written;
    }
    for f in files.values_mut() {
        f.flush()?;
    }
    fs::write(
        c.output_dir.join("relatives.json"),
        serde_json::to_vec_pretty(&summary).map_err(|e| super::cli_error(e.to_string()))?,
    )?;
    Ok(summary)
}

/// The taxon of a sample's species (`genus_species[_subspecies]`): its own,
/// else the binomial's, else the first of its genus.
fn taxon_of<'a>(taxa: &'a HashMap<String, String>, species: &str) -> Option<&'a str> {
    let words = species.split('_').collect::<Vec<_>>();
    let found = |key: &str| taxa.get(key).filter(|t| !t.is_empty());
    found(species)
        .or_else(|| words.get(..2).and_then(|w| found(&w.join("_"))))
        .or_else(|| {
            let genus = format!("{}_", words[0]);
            taxa.iter()
                .filter(|(k, t)| k.starts_with(&genus) && !t.is_empty())
                .min_by_key(|(k, _)| k.as_str())
                .map(|(_, t)| t)
        })
        .map(String::as_str)
}

/// The relatives file for a species under `rank`, as [`write_relatives`] names it.
pub fn relatives_file(dir: &Path, taxon: Option<&str>) -> PathBuf {
    dir.join(format!(
        "{}.fa",
        taxon
            .filter(|t| !t.is_empty())
            .unwrap_or(UNASSIGNED)
            .replace(['/', ' '], "_")
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_cascade_and_model_calls_per_taxon() {
        let dir = tempfile::tempdir().unwrap();
        let results = dir.path().join("results");
        let header = "family_id\toriginal_label\tdisposition\tlabel\tconfidence\tdecided_by\n";
        for (species, rows, fa) in [
            (
                "moth_a_x",
                "a1#Unknown\tUnknown\tknown_te\tDNA/hAT\tmedium\tcascade\n\
                 a2#Unknown\tUnknown\tunresolved\tLINE/L2\tmedium\tmodel\n\
                 a3#Unknown\tUnknown\tinterspersed_repeat\tLTR/Gypsy\tlow\trelative\n\
                 a4#Unknown\tUnknown\tunresolved\tUnknown\tlow\tcascade\n\
                 a5#Unknown\tUnknown\tsatellite\tSatellite\tmedium\tcascade\n\
                 a6#Unknown\tUnknown\tknown_te\tLTR/Gypsy\tmedium\tneighbour\n",
                ">a1#Unknown\nacgt\n>a2#Unknown\nAAAA\n>a3#Unknown\nCCCC\n>a4#Unknown\nGGGG\n>a5#Unknown\nTTTT\n>a6#Unknown\nGGCC\n",
            ),
            (
                "fly_b",
                "b1#LTR/Copia\tLTR/Copia\tknown_te\tLTR/Copia\thigh\tcascade\n",
                ">b1#LTR/Copia\nACGTACGT\n",
            ),
            (
                "unplaced_c",
                "c1#Unknown\tUnknown\tknown_te\tRC/Helitron\tlow\tcascade\n",
                ">c1#Unknown\nTTTT\n",
            ),
        ] {
            let d = results.join(species).join("GCA_1.1");
            fs::create_dir_all(&d).unwrap();
            fs::write(d.join("classification.tsv"), format!("{header}{rows}")).unwrap();
            fs::write(d.join("queries.fa"), fa).unwrap();
        }
        let taxa = dir.path().join("taxa.tsv");
        fs::write(
            &taxa,
            // moth_a by its subspecies, fly_b by its genus
            "species\torder\tclass\nMoth a\tLepidoptera\tInsecta\nFly c\tDiptera\tInsecta\n",
        )
        .unwrap();
        let out = dir.path().join("relatives");
        let s = write_relatives(&RelativesConfig {
            results,
            taxa: Some(taxa),
            rank: "order".into(),
            output_dir: out.clone(),
        })
        .unwrap();
        assert_eq!(s.samples, 3);
        assert_eq!(
            fs::read_to_string(out.join("Lepidoptera.fa")).unwrap(),
            ">moth_a_x:GCA_1.1|a1#DNA/hAT\nACGT\n>moth_a_x:GCA_1.1|a2#LINE/L2\nAAAA\n"
        );
        assert_eq!(
            fs::read_to_string(out.join("Diptera.fa")).unwrap(),
            ">fly_b:GCA_1.1|b1#LTR/Copia\nACGTACGT\n"
        );
        assert!(fs::read_to_string(relatives_file(&out, None))
            .unwrap()
            .starts_with(">unplaced_c:GCA_1.1|c1#RC/Helitron\n"));
        assert_eq!(s.per_taxon["Lepidoptera"], 2);
    }
}
