use super::{
    align::mafft_command,
    config::CopyRecoveryConfig,
    extract::{clip_interval, reverse_complement, FastaSequenceProvider, SequenceProvider},
    hits::{parse_blast_file, Hit},
    loci::{collapse_overlapping_loci, mark_ambiguous_strand_loci, merge_hits, rank_loci},
    manifest::{FamilyManifest, Manifest, SCHEMA_VERSION},
    search::{blastn_command, executable, makeblastdb_command},
};
use crate::{CommandRunner, Error, ErrorKind, Result};
use bio::io::fasta;
use sha2::{Digest, Sha256};
use std::{
    fs,
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
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
    fs::create_dir_all(&c.output_dir)?;
    if c.resume && !c.force {
        if let Some(existing) = load_manifest(&c.output_dir.join("manifest.json")) {
            if existing.input_fingerprint == fingerprint(c).unwrap_or_default()
                && existing.status == "success"
                && c.output_dir.join("hits.tsv").is_file()
            {
                let mut families = Vec::new();
                for entry in fs::read_dir(c.output_dir.join("families"))? {
                    let path = entry?.path();
                    if path.is_dir()
                        && path.join("manifest.json").is_file()
                        && path.join("copies.aligned.fa").is_file()
                    {
                        let id = path.file_name().unwrap().to_string_lossy().to_string();
                        families.push(FamilyRecoveryResult {
                            family_id: id,
                            directory: path,
                            n_hits: 0,
                            n_loci: 0,
                            n_selected: 0,
                        });
                    }
                }
                return Ok(CopyRecoveryResult {
                    output_dir: c.output_dir.clone(),
                    families,
                    manifest_path: c.output_dir.join("manifest.json"),
                });
            }
        }
    }
    fs::create_dir_all(c.output_dir.join("families"))?;
    let db = c.output_dir.join("genome");
    let mut cmd = makeblastdb_command(executable("makeblastdb"), &c.genome, &db);
    crate::command_runner::run_checked(r, &mut cmd, "makeblastdb")?;
    let blast = c.output_dir.join("hits.tsv");
    let mut cmd = blastn_command(executable("blastn"), &db, &c.queries, &blast, c.threads);
    crate::command_runner::run_checked(r, &mut cmd, "blastn")?;
    let hits = parse_blast_file(&blast)?;
    let provider = FastaSequenceProvider::from_file(&c.genome)?;
    let mafft_exe = executable("mafft");
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
        mark_ambiguous_strand_loci(&mut loci, c.orientation_threshold);
        rank_loci(&mut loci, c.hits_per_family, c.diversity_filter);
        let nsel = loci.iter().filter(|l| l.selected).count();
        write_hits(&dir.join("hits.tsv"), &fh)?;
        write_loci(&dir.join("loci.tsv"), &loci)?;
        let mut raw = String::new();
        let mut extended = String::new();
        let mut copies = String::from("family_id\tlocus_id\tsubject_id\tgenomic_start\tgenomic_end\tstrand\trequested_extended_start\trequested_extended_end\tactual_extended_start\tactual_extended_end\tleft_flank_actual\tright_flank_actual\treverse_complemented\tselected\textraction_status\talignment_status\texclusion_reason\n");
        let mut usable = 0usize;
        for locus in loci.iter() {
            let len = provider.contig_length(&locus.subject_id)?;
            let (rs, re) = (locus.start, locus.end);
            let (es, ee) = clip_interval(rs, re, len as usize, c.flank_length);
            let mut raw_seq = provider.sequence(&locus.subject_id, rs.min(len), re.min(len))?;
            let mut ext_seq = provider.sequence(&locus.subject_id, es as u64, ee as u64)?;
            let rc = locus.strand == super::Strand::Minus;
            if rc {
                raw_seq = reverse_complement(&raw_seq);
                ext_seq = reverse_complement(&ext_seq);
            }
            let status = if locus.selected {
                usable += 1;
                "success"
            } else {
                "not_selected"
            };
            if locus.selected {
                raw.push_str(&format!(
                    ">locus_id={} family_id={}\n{}\n",
                    locus.locus_id,
                    locus.family_id,
                    String::from_utf8_lossy(&raw_seq)
                ));
                extended.push_str(&format!(
                    ">locus_id={} family_id={}\n{}\n",
                    locus.locus_id,
                    locus.family_id,
                    String::from_utf8_lossy(&ext_seq)
                ));
            }
            copies.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                locus.family_id,
                locus.locus_id,
                locus.subject_id,
                rs,
                re,
                locus.strand.as_str(),
                es.saturating_sub(c.flank_length),
                ee.saturating_add(c.flank_length),
                es,
                ee,
                es.saturating_sub(rs as usize),
                ee.saturating_sub(re as usize),
                rc,
                locus.selected,
                status,
                if locus.selected {
                    "pending"
                } else {
                    "not_applicable"
                },
                locus.exclusion_reason.as_deref().unwrap_or("")
            ));
        }
        fs::write(dir.join("copies.tsv"), copies)?;
        fs::write(dir.join("copies.raw.fa"), raw)?;
        fs::write(dir.join("copies.extended.fa"), &extended)?;
        if usable == 0 {
            fs::write(dir.join("copies.aligned.fa"), "")?;
        } else if usable == 1 {
            fs::write(dir.join("copies.aligned.fa"), &extended)?;
        } else {
            let mut ac = mafft_command(&mafft_exe, &dir.join("copies.extended.fa"), c.threads);
            let out = crate::command_runner::run_checked(r, &mut ac, "MAFFT")?;
            fs::write(dir.join("copies.aligned.fa"), out.stdout)?;
        }
        let fm = FamilyManifest {
            schema_version: SCHEMA_VERSION.into(),
            family_id: id.clone(),
            status: if usable == 0 {
                "no_usable_copies".into()
            } else if usable == 1 {
                "single_copy".into()
            } else {
                "success".into()
            },
            raw_hits: fh.len(),
            merged_loci: loci.len(),
            selected: nsel,
            warnings: Vec::new(),
            stages: stage_map(),
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
        input_fingerprint: fingerprint(c)?,
        coordinate_convention: "BLAST 1-based inclusive; loci BED 0-based half-open".into(),
        configuration: serde_json::to_value(c)
            .map_err(|e| Error::new(ErrorKind::GenericCli(e.to_string())))?,
        warnings: Vec::new(),
        status: "success".into(),
        families: fams.len(),
        stages: stage_map(),
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
fn fingerprint(c: &CopyRecoveryConfig) -> Result<String> {
    let mut h = Sha256::new();
    for p in [&c.genome, &c.queries] {
        let b = fs::read(p)?;
        h.update(&b);
    }
    h.update(serde_json::to_vec(c).map_err(|e| Error::new(ErrorKind::GenericCli(e.to_string())))?);
    Ok(format!("{:x}", h.finalize()))
}
fn stage_map() -> std::collections::BTreeMap<String, super::StageState> {
    let mut m = std::collections::BTreeMap::new();
    for n in [
        "search", "parse", "merge", "collapse", "select", "extract", "orient", "align",
    ] {
        m.insert(
            n.into(),
            super::StageState {
                status: "complete".into(),
                started_at: None,
                completed_at: Some(
                    SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .unwrap_or_default()
                        .as_secs()
                        .to_string(),
                ),
                input_fingerprint: None,
                relevant_config: None,
                commands: Vec::new(),
                outputs: Vec::new(),
                warnings: Vec::new(),
                error: None,
            },
        );
    }
    m
}
fn load_manifest(path: &Path) -> Option<Manifest> {
    fs::read(path)
        .ok()
        .and_then(|b| serde_json::from_slice(&b).ok())
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
