use super::{
    config::CopyRecoveryConfig,
    hits::{parse_blast_file, Hit},
    loci::{collapse_overlapping_loci, merge_hits, rank_loci},
    manifest::{FamilyManifest, Manifest, SCHEMA_VERSION},
    search::{blastn_command, executable, makeblastdb_command},
};
use crate::{CommandRunner, Error, ErrorKind, Result};
use bio::io::fasta;
use std::{
    fs,
    path::{Path, PathBuf},
};
#[derive(Clone, Debug)]
pub struct FamilyRecoveryResult {
    pub family_id: String,
    pub directory: PathBuf,
    pub n_hits: usize,
    pub n_loci: usize,
    pub n_selected: usize,
}
#[derive(Clone, Debug)]
pub struct CopyRecoveryResult {
    pub output_dir: PathBuf,
    pub families: Vec<FamilyRecoveryResult>,
    pub manifest_path: PathBuf,
}
pub fn run_copy_recovery(c: &CopyRecoveryConfig) -> Result<CopyRecoveryResult> {
    run_copy_recovery_with_runner(c, &crate::RealCommandRunner)
}
pub fn run_copy_recovery_with_runner(
    c: &CopyRecoveryConfig,
    r: &dyn CommandRunner,
) -> Result<CopyRecoveryResult> {
    c.validate()
        .map_err(|e| Error::new(ErrorKind::GenericCli(e)))?;
    fs::create_dir_all(c.output_dir.join("families"))?;
    let db = c.output_dir.join("genome");
    let mut cmd = makeblastdb_command(executable("makeblastdb"), &c.genome, &db);
    crate::command_runner::run_checked(r, &mut cmd, "makeblastdb")?;
    let blast = c.output_dir.join("hits.tsv");
    let mut cmd = blastn_command(executable("blastn"), &db, &c.queries, &blast, c.threads);
    crate::command_runner::run_checked(r, &mut cmd, "blastn")?;
    let hits = parse_blast_file(&blast)?;
    let mut fams = Vec::new();
    let reader = fasta::Reader::from_file(&c.queries)?;
    for rec in reader.records() {
        let rec = rec?;
        let id = rec.id().to_string();
        let dir = c.output_dir.join("families").join(safe(&id));
        fs::create_dir_all(&dir)?;
        fs::write(
            dir.join("query.fa"),
            format!(">{}\n{}\n", id, String::from_utf8_lossy(rec.seq())),
        )?;
        let fh: Vec<Hit> = hits.iter().filter(|h| h.query_id == id).cloned().collect();
        let mut loci = collapse_overlapping_loci(merge_hits(&fh, c.max_fragment_gap as u64));
        rank_loci(&mut loci, c.hits_per_family, c.diversity_filter);
        let nsel = loci.iter().filter(|l| l.selected).count();
        write_hits(&dir.join("hits.tsv"), &fh)?;
        write_loci(&dir.join("loci.tsv"), &loci)?;
        let fm = FamilyManifest {
            schema_version: SCHEMA_VERSION.into(),
            family_id: id.clone(),
            status: "success".into(),
            raw_hits: fh.len(),
            merged_loci: loci.len(),
            selected: nsel,
            warnings: Vec::new(),
        };
        fs::write(
            dir.join("manifest.json"),
            serde_json::to_vec_pretty(&fm)
                .map_err(|e| Error::new(ErrorKind::GenericCli(e.to_string())))?,
        )?;
        fams.push(FamilyRecoveryResult {
            family_id: id,
            directory: dir,
            n_hits: fh.len(),
            n_loci: loci.len(),
            n_selected: nsel,
        });
    }
    let m = Manifest {
        schema_version: SCHEMA_VERSION.into(),
        rep_version: env!("CARGO_PKG_VERSION").into(),
        coordinate_convention: "BLAST 1-based inclusive; loci BED 0-based half-open".into(),
        configuration: serde_json::to_value(c)
            .map_err(|e| Error::new(ErrorKind::GenericCli(e.to_string())))?,
        warnings: Vec::new(),
        status: "success".into(),
        families: fams.len(),
    };
    let mp = c.output_dir.join("manifest.json");
    fs::write(
        &mp,
        serde_json::to_vec_pretty(&m)
            .map_err(|e| Error::new(ErrorKind::GenericCli(e.to_string())))?,
    )?;
    Ok(CopyRecoveryResult {
        output_dir: c.output_dir.clone(),
        families: fams,
        manifest_path: mp,
    })
}
fn safe(s: &str) -> String {
    let x: String = s
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' || c == '-' {
                c
            } else {
                '_'
            }
        })
        .collect();
    if x.is_empty() {
        "family".into()
    } else {
        x
    }
}
fn write_hits(p: &Path, h: &[Hit]) -> Result<()> {
    let mut s: String="family_id\tquery_id\tsubject_id\tquery_start\tquery_end\tquery_length\tsubject_start\tsubject_end\tstrand\talignment_length\tmatches\tmismatches\tgap_opens\tpercent_identity\tquery_coverage\tbitscore\tevalue\n".to_string();
    for x in h {
        s.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            x.family_id,
            x.query_id,
            x.subject_id,
            x.query_start,
            x.query_end,
            x.query_length,
            x.subject_start,
            x.subject_end,
            x.strand.as_str(),
            x.alignment_length,
            x.matches,
            x.mismatches,
            x.gap_opens,
            x.percent_identity,
            x.query_coverage,
            x.bitscore,
            x.evalue
        ));
    }
    fs::write(p, s)?;
    Ok(())
}
fn write_loci(p: &Path, l: &[super::Locus]) -> Result<()> {
    let mut s: String="family_id\tlocus_id\tsubject_id\tstart\tend\tstrand\tn_fragments\tbest_bitscore\tcombined_query_coverage\trepresentative_score\tselected\texclusion_reason\n".into();
    for x in l {
        s.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            x.family_id,
            x.locus_id,
            x.subject_id,
            x.start,
            x.end,
            x.strand.as_str(),
            x.n_fragments,
            x.best_bitscore,
            x.combined_query_coverage,
            x.representative_score,
            x.selected,
            x.exclusion_reason.as_deref().unwrap_or("")
        ));
    }
    fs::write(p, s)?;
    Ok(())
}
