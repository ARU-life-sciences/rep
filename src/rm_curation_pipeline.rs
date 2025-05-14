// give proper acknowledgement to the original source code and provide a link to the original source code.
use std::{
    cmp::{max, min},
    collections::{HashMap, HashSet},
    fs::{self, OpenOptions},
    io::BufWriter,
    path::{Path, PathBuf},
    process::Command,
};

use crate::{parse_blast::BlastRecord, CliArgs, Error, ErrorKind, Result, DATA, INTERMEDIATE};

use anyhow::Context;
use bio::io::fasta::{self, Writer};

// TODO: before we start, have to make sure there
// are no special characters in the fasta headers
// `/`. Do this here.

pub fn rm_curation_pipeline(matches: CliArgs) -> Result<()> {
    let mut data_path = matches.curation_only.clone().unwrap();
    data_path.push(DATA);

    // TODO: check this is correct. The aim here is to use the
    // fasta we copied to the data directory. So we need its name.
    let cli_matches_fasta_file = matches.fasta_file.clone();
    // get the basename and push this onto the data path

    let fasta_name = match cli_matches_fasta_file.file_name() {
        Some(b) => {
            // if the fasta name is gzipped, we remove the .gz extension
            // as this is already dealt with in previous steps

            if b.to_string_lossy().ends_with(".gz") {
                PathBuf::from(b.to_string_lossy().split(".gz").collect::<Vec<_>>()[0])
            } else {
                PathBuf::from(b)
            }
        }
        None => {
            return Err(Error::new(ErrorKind::GenericCli(
                "Could not get fasta name".into(),
            )))
        }
    };

    eprintln!("The fasta file is {:?}", fasta_name);

    // the original data path
    let data_path_clone = data_path.clone();
    // the genome fasta name
    data_path.push(fasta_name.clone());

    // first of all we need to run "makeblastdb"
    // with the input fasta,
    // but we want the dbtype to be nucleotide and
    // parse the seqids
    eprintln!("Running makeblastdb in dir: {}", data_path_clone.display());
    let makeblastdb =
        std::process::Command::new("/software/team301/ncbi-blast-2.16.0+/bin/makeblastdb")
            .current_dir(data_path_clone)
            .args(["-in", fasta_name.to_str().unwrap()])
            // this is what is written here:
            // https://github.com/ValentinaPeona/TardigraTE/blob/main/Practicals/Practical3.md
            // but we can modify the name later
            .args(["-out", fasta_name.to_str().unwrap()])
            .args(["-dbtype", "nucl"])
            .arg("-parse_seqids")
            .spawn()
            .with_context(|| "makeblastdb command failed".to_string())?;

    let _blastdb = makeblastdb.wait_with_output()?;

    // this path comes from the `matches` input
    let rmdl_library = matches.curation_rmdl_library.clone().unwrap();

    eprintln!("Running the RM curation pipeline");
    run_rmdl_curation_pipeline(matches, fasta_name, rmdl_library)?;

    Ok(())
}

// the actual pipeline
fn run_rmdl_curation_pipeline(
    matches: CliArgs,
    genome_fasta_name: PathBuf,
    rmdl_library: PathBuf,
) -> Result<()> {
    // FIXME: this path is used a bunch of times, maybe we can
    // factor it out.
    // the path to the genome, and also the blast database
    let mut genome_and_blastdb_path = matches.curation_only.clone().unwrap();
    genome_and_blastdb_path.push(DATA);
    genome_and_blastdb_path.push(genome_fasta_name.clone());

    // # Make folders blastn and aligned
    // we want these folders to be in the intermediate directory

    let blastn_dir = matches
        .curation_only
        .clone()
        .unwrap()
        .join(INTERMEDIATE)
        .join("blastn");
    fs::create_dir_all(blastn_dir.clone())?;
    let aligned_dir = matches
        .curation_only
        .clone()
        .unwrap()
        .join(INTERMEDIATE)
        .join("aligned");

    fs::create_dir_all(aligned_dir.clone())?;

    let flank = 2000;
    let maxhitdist = 10000;
    let minfrac = 0.8;
    // let digits = 5;
    let hits = 20;

    // TODO: keep in mind these will probably need to
    // be in context of the matches.configure
    let temp_blast_out = matches
        .curation_only
        .clone()
        .unwrap()
        .join(INTERMEDIATE)
        .join("tempBlastOut.txt");

    let temp_map_names = matches
        .curation_only
        .clone()
        .unwrap()
        .join(INTERMEDIATE)
        .join("tempMapNames.txt");

    eprintln!("The blast output file is {:?}", temp_blast_out);
    eprintln!("The map names file is {:?}", temp_map_names);

    // the blastdb is the fasta file we created
    // without the extension
    let blastdb = genome_fasta_name.clone();
    // 1. blast all files in the repeatmasked directory
    blast_repeatmasked(
        matches.clone(),
        rmdl_library.clone(),
        blastdb,
        temp_blast_out.clone(),
        hits,
    )?;

    // 2. Find hits from assembly (genome) and add original query
    find_hits_from_assembly(
        matches.clone(),
        blastn_dir.clone(),
        maxhitdist,
        minfrac,
        genome_fasta_name,
        temp_blast_out,
        flank,
    )?;

    // 3. Run alignments with MAFFT
    let mafft_log = matches
        .curation_only
        .clone()
        .unwrap()
        .join(INTERMEDIATE)
        .join("tempMafft.txt");

    run_mafft_alignments(&blastn_dir, &aligned_dir, &mafft_log)?;

    Ok(())
}

// 1. blast all files in the repeatmasked directory
fn blast_repeatmasked(
    matches: CliArgs,
    rmdl_library: PathBuf,
    blastdb: PathBuf,
    temp_blast: PathBuf,
    _hits: i32,
) -> Result<()> {
    let mut data_path = matches.curation_only.clone().unwrap();
    data_path.push(DATA);

    // create the blastdb and the rmdl_library variables
    // these are relative to the data path
    let blastdb_data_path = data_path.clone();
    let blastdb = blastdb_data_path.join(blastdb); // name of the blast database
    let rmdl_library_data_path = data_path.clone();
    let rmdl_library = rmdl_library_data_path.join(rmdl_library); // name of the query file

    eprintln!(
        "rep :: blast_repeatmasked :: BLAST all sequences in {:?} against {:?} ",
        rmdl_library, blastdb
    );

    let blastn_command = std::process::Command::new(
        "/software/team301/ncbi-blast-2.16.0+/bin/blastn",
    )
    .args(["-db", &blastdb.to_string_lossy()])
    .args(["-query", &rmdl_library.to_string_lossy()])
    // 7 includes header lines
    .args([
        "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
    ])
    .arg("-evalue")
    .arg("10e-10")
    .args(["-out", &temp_blast.to_string_lossy()])
    .status()?;

    if !blastn_command.success() {
        return Err(Error::new(ErrorKind::GenericCli("BLASTN failed".into())));
    }

    eprintln!("BLAST complete");

    Ok(())
}

// 2. Find hits from assembly (genome) and add original query
// &step2($FASTA, $blastdir, $maxhitdist, $minfrac, $ASSEMBLY);

fn find_hits_from_assembly(
    matches: CliArgs,
    blastn_dir: PathBuf,
    maxhitdist: i32,
    minfrac: f64,
    genome_fasta_name: PathBuf,
    temp_blast_out: PathBuf,
    flank: u64,
) -> Result<()> {
    // Construct genome FASTA path
    let mut genome_path = matches.curation_only.clone().unwrap();
    genome_path.push(DATA);
    genome_path.push(genome_fasta_name);

    // Load genome into memory once
    let genome = load_genome(&genome_path)?;

    // Clear existing output
    for entry in fs::read_dir(blastn_dir.clone())? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().map(|e| e == "fa").unwrap_or(false)
            || path
                .components()
                .any(|c| c.as_os_str().to_string_lossy().contains("txt"))
        {
            fs::remove_file(&path)?;
        }
    }

    eprintln!("rep :: find_hits_from_assembly :: Processing BLAST hits...");

    // Load and group BLAST hits
    let blast_hits = BlastRecord::from_file(temp_blast_out)?;
    for h in &blast_hits.0 {
        eprintln!(
            "HIT: qseqid={}, sseqid={}, evalue={}, sstart={}, send={}",
            h.qseqid, h.sseqid, h.evalue, h.sstart, h.send
        );
    }
    let repeat_ids: HashSet<_> = blast_hits.0.iter().map(|h| h.qseqid.clone()).collect();

    for repeat_id in repeat_ids {
        let outfile = blastn_dir.join(format!("{}.fa", sanitize_filename(&repeat_id)));

        let mut filtered_hits = blast_hits.filter_by_query_name(&repeat_id);
        filtered_hits.sort_by_evalue();
        eprintln!(
            "Processing query {} with {} hits",
            repeat_id,
            filtered_hits.0.len()
        );
        let top_hits = filtered_hits.top_n(20).filter_unique_combinations();

        for hit in top_hits.0 {
            let genome_id = &hit.sseqid;

            let mut strand_plus = 0;
            let mut strand_minus = 0;
            let mut start = 0;
            let mut stop = 0;

            let mut per_pair_hits = blast_hits.filter_by_query_subject(&repeat_id, genome_id);
            per_pair_hits.sort_by_alignment_positions();

            for (i, h) in per_pair_hits.0.iter().enumerate() {
                let hs = h.sstart;
                let he = h.send;

                if i == 0 {
                    start = min(hs, he);
                    stop = max(hs, he);
                    if he > hs {
                        strand_plus += 1;
                    } else {
                        strand_minus += 1;
                    }
                } else if min(hs, he) - stop < maxhitdist as u64 {
                    start = min(start, min(hs, he));
                    stop = max(stop, max(hs, he));
                    if he > hs {
                        strand_plus += 1;
                    } else {
                        strand_minus += 1;
                    }
                } else {
                    let direction = if strand_plus >= strand_minus {
                        '+'
                    } else {
                        '-'
                    };
                    if (max(strand_plus, strand_minus) as f64 / (strand_plus + strand_minus) as f64)
                        < minfrac
                    {
                        eprintln!("Ambiguous orientation for {} on {}", repeat_id, genome_id);
                    }

                    extract_seq(
                        &genome,
                        genome_id,
                        start.saturating_sub(flank),
                        stop + flank,
                        direction,
                        outfile.clone(),
                    )
                    .with_context(|| {
                        format!(
                            "Failed to extract {} from {}:{}-{}",
                            repeat_id, genome_id, start, stop
                        )
                    })?;

                    // reset for next cluster
                    start = min(hs, he);
                    stop = max(hs, he);
                    if he > hs {
                        strand_plus = 1;
                        strand_minus = 0;
                    } else {
                        strand_plus = 0;
                        strand_minus = 1;
                    }
                }
            }

            // Final segment
            let direction = if strand_plus >= strand_minus {
                '+'
            } else {
                '-'
            };
            if (max(strand_plus, strand_minus) as f64 / (strand_plus + strand_minus) as f64)
                < minfrac
            {
                eprintln!("Ambiguous orientation for {} on {}", repeat_id, genome_id);
            }

            extract_seq(
                &genome,
                genome_id,
                start.saturating_sub(flank),
                stop + flank,
                direction,
                outfile.clone(),
            )
            .with_context(|| {
                format!(
                    "Failed to extract {} from {}:{}-{}",
                    repeat_id, genome_id, start, stop
                )
            })?;
        }
    }

    Ok(())
}

fn run_mafft_alignments(blast_dir: &Path, align_dir: &Path, log_path: &Path) -> Result<()> {
    fs::create_dir_all(align_dir)?;

    let log_file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(log_path)
        .with_context(|| format!("Failed to open file {:?}", log_path))?;

    for entry in fs::read_dir(blast_dir)? {
        let entry = entry?;
        let path = entry.path();

        if path.extension().map(|e| e == "fa").unwrap_or(false) {
            if fs::metadata(&path)?.len() == 0 {
                eprintln!("Skipping empty file: {:?}", path);
                continue;
            }

            if let Some(file_name) = path.file_name() {
                let output = align_dir.join(file_name);

                eprintln!("Running MAFFT on {:?}", file_name);

                let status =
                    Command::new("/software/team301/mafft-7.525-with-extensions/core/mafft")
                        .arg("--ep")
                        .arg("0.0")
                        .arg("--genafpair")
                        .arg("--maxiterate")
                        .arg("1000")
                        .arg("--thread")
                        .arg("10")
                        .arg("--adjustdirection")
                        .arg(path.to_str().unwrap())
                        .stdout(fs::File::create(&output)?)
                        .stderr(log_file.try_clone()?)
                        .status()?;

                if !status.success() {
                    return Err(Error::new(ErrorKind::GenericCli(format!(
                        "MAFFT failed for {:?}",
                        path
                    ))));
                }
            }
        }
    }

    Ok(())
}

pub fn load_genome(fasta: &PathBuf) -> Result<HashMap<String, Vec<u8>>> {
    let mut genome = HashMap::new();
    let reader = fasta::Reader::from_file(fasta)?;
    for rec in reader.records() {
        let rec = rec?;
        genome.insert(rec.id().to_string(), rec.seq().to_vec());
    }
    Ok(genome)
}

// use an iterator approach to extract sequences
fn extract_seq(
    genome: &HashMap<String, Vec<u8>>,
    query: &str,
    start: u64,
    stop: u64,
    direction: char,
    outfile: PathBuf,
) -> Result<()> {
    let file = OpenOptions::new()
        .append(true)
        .create(true)
        .open(outfile.clone())
        .with_context(|| format!("Failed to open file {:?}", outfile))?;
    let mut writer = Writer::new(BufWriter::new(file));

    if let Some(seq) = genome.get(query) {
        let start = start.saturating_sub(1);
        let end = stop.min(seq.len() as u64);
        let mut substr = seq[start as usize..end as usize].to_vec();
        if direction == '-' {
            substr = bio::alphabets::dna::revcomp(substr);
        }
        writer.write(query, None, &substr)?;
    } else {
        eprintln!("WARN: sequence {} not found in genome", query);
    }

    Ok(())
}

fn sanitize_filename(s: &str) -> String {
    s.replace(['/', '\\', ' ', '#', ':'], "_")
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;

    #[test]
    fn test_load_genome() {
        let path = PathBuf::from("./test/data/genome.fa");
        let genome = load_genome(&path).unwrap();
        // one scaffold
        assert_eq!(genome.len(), 1);
        assert_eq!(
            genome["chr1"],
            b"ATGCGTACGTAGCTAGCTGACTGATCGATCGTAGCTAGCTAGCTGATCGTACGTAGCTAG"
        );
    }

    #[test]
    fn test_extract_seq() {
        let genome = load_genome(&PathBuf::from("./test/data/genome.fa")).unwrap();
        let output = PathBuf::from("./test/test_extract_seq.fa");

        extract_seq(&genome, "scaffold1", 1, 16, '+', output.clone()).unwrap();

        let contents = std::fs::read_to_string(output).unwrap();
        assert!(contents.contains(">scaffold1"));
        assert!(contents.contains("ACGTACGTACGTACGT"));
    }

    #[test]
    fn test_full_pipeline() -> Result<()> {
        let genome = PathBuf::from("./test/data/genome2.fa");
        let repeat = PathBuf::from("rmdl2.fa");

        let args = CliArgs {
            fasta_file: genome,
            configure: None,
            database: None,
            rmo_threads: 1,
            rma_threads: 1,
            rma_only: false,
            curation_only: Some(PathBuf::from("./test")),
            curation_rmdl_library: Some(repeat),
            verbose: false,
        };

        rm_curation_pipeline(args)?; // your entrypoint

        Ok(())
    }
}
