use std::{
    collections::{BTreeMap, HashSet},
    fs::{self, File},
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

use bio::{alphabets::dna, io::fasta};

use crate::{
    command_runner::run_checked,
    configured_executable,
    parse_blast::{BlastRecord, BlastTable},
    repeatmasker::find_consensus_library,
    CliArgs, CommandRunner, Error, ErrorKind, Result, BLASTN_ENV, DATA, MAFFT_ENV, MAKEBLASTDB_ENV,
};

const CURATION: &str = "curation";
const BLAST_HITS: &str = "blastn";
const ALIGNMENTS: &str = "aligned";
const BLAST_OUTPUT: &str = "blast_hits.tsv";
const DEFAULT_FLANK: u64 = 2_000;
const DEFAULT_MAX_HIT_DISTANCE: u64 = 10_000;
const DEFAULT_MIN_ORIENTATION_FRACTION: f64 = 0.8;
const DEFAULT_HITS_PER_QUERY: usize = 20;

#[derive(Clone, Copy, Debug)]
pub struct CurationConfig {
    pub flank: u64,
    pub max_hit_distance: u64,
    pub min_orientation_fraction: f64,
    pub hits_per_query: usize,
}

impl Default for CurationConfig {
    fn default() -> Self {
        Self {
            flank: DEFAULT_FLANK,
            max_hit_distance: DEFAULT_MAX_HIT_DISTANCE,
            min_orientation_fraction: DEFAULT_MIN_ORIENTATION_FRACTION,
            hits_per_query: DEFAULT_HITS_PER_QUERY,
        }
    }
}

#[derive(Debug)]
struct HitRegion {
    start: u64,
    end: u64,
    plus: usize,
    minus: usize,
}

pub fn run_curation_via_copy_recovery(matches: CliArgs, runner: &dyn CommandRunner) -> Result<()> {
    let configure = matches.configure_path()?;
    let queries = find_consensus_library(&configure.join(DATA).join("repeatmodeler"))?;
    let genome = matches
        .curation_target
        .clone()
        .unwrap_or_else(|| matches.fasta_file.clone());
    let out = configure.join(CURATION).join("recovery");
    let cfg = crate::CopyRecoveryConfig {
        genome,
        queries,
        output_dir: out.clone(),
        flank_length: 2000,
        hits_per_family: 20,
        max_fragment_gap: 10000,
        orientation_threshold: 0.8,
        threads: matches.rmo_threads as usize,
        search_backend: crate::SearchBackend::Blastn,
        blast_task: crate::copy_recovery::default_blast_task(),
        evalue: crate::copy_recovery::default_evalue(),
        aligner: crate::Aligner::Mafft,
        force: false,
        resume: false,
        diversity_filter: false,
        reputils: None,
    };
    let result = crate::copy_recovery::run_copy_recovery_with_runner(&cfg, runner)?;
    let legacy_hits = configure.join(CURATION).join(BLAST_HITS);
    let legacy_align = configure.join(CURATION).join(ALIGNMENTS);
    fs::create_dir_all(&legacy_hits)?;
    fs::create_dir_all(&legacy_align)?;
    for family in result.families {
        let name = format!("{}.fa", safe_file_name(&family.family_id));
        fs::copy(
            family.directory.join("copies.extended.fa"),
            legacy_hits.join(&name),
        )?;
        fs::copy(
            family.directory.join("copies.aligned.fa"),
            legacy_align.join(&name),
        )?;
    }
    Ok(())
}

pub fn run_curation(matches: CliArgs, runner: &dyn CommandRunner) -> Result<()> {
    run_curation_with_config(matches, runner, CurationConfig::default())
}

pub fn run_curation_with_config(
    matches: CliArgs,
    runner: &dyn CommandRunner,
    parameters: CurationConfig,
) -> Result<()> {
    if parameters.hits_per_query == 0 {
        return Err(Error::new(ErrorKind::GenericCli(
            "Curation hits_per_query must be greater than zero".into(),
        )));
    }
    if !(0.0..=1.0).contains(&parameters.min_orientation_fraction) {
        return Err(Error::new(ErrorKind::GenericCli(
            "Curation min_orientation_fraction must be between 0 and 1".into(),
        )));
    }

    let configure = matches.configure_path()?;
    let curation_path = configure.join(CURATION);
    let blast_hits_path = curation_path.join(BLAST_HITS);
    let alignment_path = curation_path.join(ALIGNMENTS);
    fs::create_dir_all(&blast_hits_path)?;
    fs::create_dir_all(&alignment_path)?;

    let target_source = matches
        .curation_target
        .as_deref()
        .unwrap_or(&matches.fasta_file);
    let target_fasta = prepare_target_fasta(target_source, &curation_path)?;

    let repeatmodeler_path = configure.join(DATA).join("repeatmodeler");
    let consensus_path = find_consensus_library(&repeatmodeler_path)?;
    let database_path = curation_path.join("genome");

    let mut makeblastdb = Command::new(configured_executable(MAKEBLASTDB_ENV, "makeblastdb"));
    makeblastdb
        .arg("-in")
        .arg(&target_fasta)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&database_path)
        .arg("-parse_seqids");
    run_checked(runner, &mut makeblastdb, "makeblastdb")?;

    let blast_output_path = curation_path.join(BLAST_OUTPUT);
    let mut blastn = Command::new(configured_executable(BLASTN_ENV, "blastn"));
    blastn
        .arg("-db")
        .arg(&database_path)
        .arg("-query")
        .arg(&consensus_path)
        .arg("-outfmt")
        .arg("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore")
        .arg("-evalue")
        .arg("1e-9")
        .arg("-max_target_seqs")
        .arg(parameters.hits_per_query.to_string())
        .arg("-out")
        .arg(&blast_output_path);
    run_checked(runner, &mut blastn, "blastn")?;

    let blast_table = BlastRecord::from_file(&blast_output_path)?;
    let subjects = blast_table
        .records()
        .iter()
        .flat_map(|record| {
            [
                record.sseqid.clone(),
                normalized_sequence_id(&record.sseqid).to_owned(),
            ]
        })
        .collect::<HashSet<_>>();
    let genome = read_fasta(&target_fasta, Some(&subjects))?;
    let consensus = read_fasta(&consensus_path, None)?;

    let mut output_names = HashSet::new();
    for (query_id, query_sequence) in consensus {
        let query_hits = blast_table.filter_by_query_name(&query_id);
        let copies = extract_related_copies(&query_id, query_hits, &genome, parameters)?;
        let file_name = format!("{}.fa", safe_file_name(&query_id));
        if !output_names.insert(file_name.clone()) {
            return Err(Error::new(ErrorKind::GenericCli(format!(
                "Multiple consensus identifiers map to output file {file_name}"
            ))));
        }
        let unaligned_path = blast_hits_path.join(&file_name);
        write_query_and_copies(&unaligned_path, &query_id, &query_sequence, &copies)?;

        let aligned_path = alignment_path.join(file_name);
        let mut mafft = Command::new(configured_executable(MAFFT_ENV, "mafft"));
        mafft
            .arg("--ep")
            .arg("0.0")
            .arg("--genafpair")
            .arg("--maxiterate")
            .arg("1000")
            .arg("--thread")
            .arg(matches.rmo_threads.to_string())
            .arg("--adjustdirection")
            .arg(&unaligned_path);
        let output = run_checked(runner, &mut mafft, "MAFFT")?;
        fs::write(aligned_path, output.stdout)?;
    }

    Ok(())
}

fn prepare_target_fasta(source: &Path, curation_path: &Path) -> Result<PathBuf> {
    if source.extension().and_then(|extension| extension.to_str()) != Some("gz") {
        return Ok(source.to_path_buf());
    }

    let target = curation_path.join("curation_target.fa");
    let mut child = Command::new("gunzip")
        .arg("-c")
        .arg(source)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;
    let mut stdout = child.stdout.take().ok_or_else(|| {
        Error::new(ErrorKind::GenericCli(
            "Could not capture gunzip output for curation target".into(),
        ))
    })?;
    let mut output_file = File::create(&target)?;
    std::io::copy(&mut stdout, &mut output_file)?;
    let output = child.wait_with_output()?;
    if !output.status.success() {
        let _ = fs::remove_file(&target);
        return Err(Error::new(ErrorKind::GenericCli(format!(
            "Could not decompress curation target: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        ))));
    }
    Ok(target)
}

fn read_fasta(
    path: &Path,
    selected_ids: Option<&HashSet<String>>,
) -> Result<BTreeMap<String, Vec<u8>>> {
    let reader = fasta::Reader::from_file(path)?;
    let mut sequences = BTreeMap::new();
    let mut seen = HashSet::new();
    for record in reader.records() {
        let record = record?;
        if !seen.insert(record.id().to_owned()) {
            return Err(Error::new(ErrorKind::GenericCli(format!(
                "Duplicate FASTA identifier {} in {}",
                record.id(),
                path.display()
            ))));
        }
        if selected_ids.is_none_or(|ids| ids.contains(record.id())) {
            sequences.insert(record.id().to_owned(), record.seq().to_vec());
        }
    }
    Ok(sequences)
}

fn normalized_sequence_id(id: &str) -> &str {
    let mut fields = id.split('|');
    match (fields.next(), fields.next()) {
        (Some("emb" | "gb" | "dbj" | "ref" | "lcl"), Some(accession)) if !accession.is_empty() => {
            accession
        }
        _ => id,
    }
}

fn extract_related_copies(
    query_id: &str,
    mut hits: BlastTable,
    genome: &BTreeMap<String, Vec<u8>>,
    parameters: CurationConfig,
) -> Result<Vec<(String, Vec<u8>)>> {
    hits.sort_by_evalue();
    let hits = hits.top_n(parameters.hits_per_query);
    let mut by_subject: BTreeMap<String, Vec<BlastRecord>> = BTreeMap::new();
    for hit in hits.into_records() {
        by_subject.entry(hit.sseqid.clone()).or_default().push(hit);
    }

    let mut copies = Vec::new();
    for (subject, subject_hits) in by_subject {
        let sequence = genome
            .get(&subject)
            .or_else(|| genome.get(normalized_sequence_id(&subject)))
            .ok_or_else(|| {
                Error::new(ErrorKind::GenericCli(format!(
                    "BLAST subject {subject} for query {query_id} is absent from the genome FASTA"
                )))
            })?;

        for region in merge_hits(subject_hits, parameters.max_hit_distance) {
            let total = region.plus + region.minus;
            let reverse = region.minus > region.plus;
            let major_fraction = region.plus.max(region.minus) as f64 / total as f64;
            if major_fraction < parameters.min_orientation_fraction {
                eprintln!(
                    "Ambiguous orientation for {query_id} hit on {subject}: {} +, {} -",
                    region.plus, region.minus
                );
            }

            let sequence_len = sequence.len() as u64;
            let start = region.start.saturating_sub(parameters.flank).max(1);
            let end = region
                .end
                .saturating_add(parameters.flank)
                .min(sequence_len);
            if start > end || end == 0 {
                return Err(Error::new(ErrorKind::GenericCli(format!(
                    "BLAST coordinates {start}-{end} are outside subject {subject}"
                ))));
            }

            let mut extracted = sequence[(start - 1) as usize..end as usize].to_vec();
            let direction = if reverse {
                extracted = dna::revcomp(&extracted);
                "minus"
            } else {
                "plus"
            };
            let output_subject = normalized_sequence_id(&subject);
            copies.push((
                format!("{output_subject}_{start}-{end}_{direction}"),
                extracted,
            ));
        }
    }
    Ok(copies)
}

fn merge_hits(mut hits: Vec<BlastRecord>, max_distance: u64) -> Vec<HitRegion> {
    hits.sort_by_key(|hit| hit.sstart.min(hit.send));
    let mut regions: Vec<HitRegion> = Vec::new();

    for hit in hits {
        let start = hit.sstart.min(hit.send);
        let end = hit.sstart.max(hit.send);
        let is_plus = hit.send >= hit.sstart;

        if let Some(region) = regions.last_mut() {
            if start.saturating_sub(region.end) <= max_distance {
                region.start = region.start.min(start);
                region.end = region.end.max(end);
                region.plus += usize::from(is_plus);
                region.minus += usize::from(!is_plus);
                continue;
            }
        }

        regions.push(HitRegion {
            start,
            end,
            plus: usize::from(is_plus),
            minus: usize::from(!is_plus),
        });
    }
    regions
}

fn write_query_and_copies(
    path: &Path,
    query_id: &str,
    query_sequence: &[u8],
    copies: &[(String, Vec<u8>)],
) -> Result<()> {
    let mut writer = fasta::Writer::to_file(path)?;
    writer.write_record(&fasta::Record::with_attrs(query_id, None, query_sequence))?;
    for (id, sequence) in copies {
        writer.write_record(&fasta::Record::with_attrs(id, None, sequence))?;
    }
    Ok(())
}

fn safe_file_name(id: &str) -> String {
    id.chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() || matches!(character as u32, 45 | 46 | 95) {
                character
            } else {
                char::from(95_u8)
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(start: u64, end: u64) -> BlastRecord {
        BlastRecord {
            qseqid: "query".into(),
            sseqid: "subject".into(),
            pident: 100.0,
            length: start.abs_diff(end) + 1,
            mismatch: 0,
            gapopen: 0,
            qstart: 1,
            qend: 10,
            sstart: start,
            send: end,
            evalue: 1e-20,
            bitscore: 100.0,
        }
    }

    #[test]
    fn merges_nearby_hits_and_tracks_orientation() {
        let regions = merge_hits(vec![hit(100, 150), hit(225, 175), hit(500, 550)], 50);
        assert_eq!(regions.len(), 2);
        assert_eq!((regions[0].start, regions[0].end), (100, 225));
        assert_eq!((regions[0].plus, regions[0].minus), (1, 1));
    }

    #[test]
    fn normalizes_blast_database_sequence_ids() {
        assert_eq!(
            normalized_sequence_id("emb|CBAADT010000001.1|"),
            "CBAADT010000001.1"
        );
        assert_eq!(normalized_sequence_id("plain_id"), "plain_id");
    }

    #[test]
    fn sanitizes_query_ids_for_output_paths() {
        assert_eq!(safe_file_name("family/1 (DNA)"), "family_1__DNA_");
    }
    struct MockRunner;

    impl CommandRunner for MockRunner {
        fn run(&self, cmd: &mut Command) -> Result<std::process::Output> {
            use std::os::unix::process::ExitStatusExt;

            let program = cmd.get_program().to_string_lossy();
            let args = cmd
                .get_args()
                .map(|argument| argument.to_string_lossy().to_string())
                .collect::<Vec<_>>();

            if program.ends_with("blastn") {
                let output_index = args.iter().position(|argument| argument == "-out").unwrap();
                fs::write(
                    &args[output_index + 1],
                    "family1\temb|chr1|\t100.0\t4\t0\t0\t1\t4\t5\t8\t1e-20\t50.0\n",
                )?;
            }

            let stdout = if program.ends_with("mafft") {
                b">family1\nACGT\n>chr1_3-10_plus\nAACCGGTT\n".to_vec()
            } else {
                Vec::new()
            };
            Ok(std::process::Output {
                status: std::process::ExitStatus::from_raw(0),
                stdout,
                stderr: Vec::new(),
            })
        }
    }

    #[test]
    fn runs_blast_extraction_and_alignment_workflow() {
        let temp = tempfile::tempdir().unwrap();
        let genome = temp.path().join("genome.fa");
        fs::write(&genome, ">chr1\nAAAACCCCGGGGTTTT\n").unwrap();

        let repeatmodeler = temp.path().join(DATA).join("repeatmodeler").join("RM_1");
        fs::create_dir_all(&repeatmodeler).unwrap();
        fs::write(
            repeatmodeler.join("consensi.fa.classified"),
            ">family1 (DNA)\nACGT\n",
        )
        .unwrap();

        let args = CliArgs {
            fasta_file: genome,
            curation_target: None,
            configure: Some(temp.path().to_path_buf()),
            database: Some("genome".into()),
            rmo_threads: 1,
            rma_threads: 1,
            rma_only: false,
            curation_only: false,
            verbose: false,
            ltr_struct: false,
            refdb: None,
            pfam: None,
        };
        let config = CurationConfig {
            flank: 2,
            max_hit_distance: 10,
            min_orientation_fraction: 0.8,
            hits_per_query: 20,
        };

        run_curation_with_config(args, &MockRunner, config).unwrap();

        let extracted = fs::read_to_string(
            temp.path()
                .join(CURATION)
                .join(BLAST_HITS)
                .join("family1.fa"),
        )
        .unwrap();
        assert!(extracted.contains(">family1"));
        assert!(extracted.contains(">chr1_3-10_plus"));
        assert!(temp
            .path()
            .join(CURATION)
            .join(ALIGNMENTS)
            .join("family1.fa")
            .is_file());
    }
}
