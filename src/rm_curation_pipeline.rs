// give proper acknowledgement to the original source code and provide a link to the original source code.
use std::{
    cmp::{max, min},
    fs::{self, OpenOptions},
    io::{self, BufWriter},
    path::PathBuf,
};

use crate::{parse_blast::BlastRecord, CliArgs, Error, ErrorKind, Result, DATA, INTERMEDIATE};

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
    eprintln!("Running makeblastdb");
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
            .spawn()?;

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
    fs::create_dir_all(aligned_dir)?;

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

    // should we be making the mafft here?
    // let temp_mafft = matches
    //     .curation_only
    //     .clone()
    //     .unwrap()
    //     .join(INTERMEDIATE)
    //     .join("tempMafft.txt");

    // the blastdb is the fasta file we created
    // without the extension
    let blastdb = genome_fasta_name.clone();
    // 1. blast all files in the repeatmasked directory
    blast_repeatmasked(
        matches.clone(),
        rmdl_library.clone(),
        blastdb,
        temp_blast_out.clone(),
        temp_map_names,
        hits,
    )?;

    // 2. Find hits from assembly (genome) and add original query
    find_hits_from_assembly(
        matches,
        blastn_dir,
        maxhitdist,
        minfrac,
        genome_fasta_name,
        temp_blast_out,
        flank,
    )?;

    Ok(())
}

// 1. blast all files in the repeatmasked directory
fn blast_repeatmasked(
    matches: CliArgs,
    rmdl_library: PathBuf,
    blastdb: PathBuf,
    temp_blast: PathBuf,
    temp_map_names: PathBuf,
    // FIXME: do I need this??
    _hits: i32,
) -> Result<()> {
    let mut data_path = matches.curation_only.clone().unwrap();
    data_path.push(DATA);

    // remove the old temp map names file if it's present
    // FIXME: what is this?
    if temp_map_names.exists() {
        fs::remove_file(&temp_map_names)?;
    }

    // create a file with the name of the temp_map_names
    // FIXME: what's the point of this??
    // let temp_map_names_file = fs::File::create(&temp_map_names)?;

    // TODO: check this is okay not to use `current_dir()`
    // I think as all inputs and outputs have paths specified, this is
    // okay.
    //
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

    let blastn_command =
        std::process::Command::new("/software/team301/ncbi-blast-2.16.0+/bin/blastn")
            .args(["-db", &blastdb.to_string_lossy()])
            .args(["-query", &rmdl_library.to_string_lossy()])
            // 7 includes header lines
            .args(["-outfmt", "7"])
            .arg("-evalue")
            .arg("10e-10")
            .args(["-out", &temp_blast.to_string_lossy()])
            .spawn()?;

    eprintln!("blastn command succeeded");

    let _ = blastn_command.wait_with_output()?;

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
    // the path to the genome, and also the blast database
    let mut genome_path = matches.curation_only.clone().unwrap();
    genome_path.push(DATA);
    genome_path.push(genome_fasta_name);

    // Remove any already existing fasta files in the blast dir
    for entry in fs::read_dir(blastn_dir.clone())? {
        let entry = entry?;
        let path = entry.path();
        // assuming all fastas end with "fa"
        if path.extension().map(|e| e == "fa").unwrap_or(false)
            || path
                .components()
                // this should kill any files with "txt" in them
                .any(|c| c.as_os_str().to_string_lossy().contains("txt"))
        {
            fs::remove_file(&path)?;
        }
    }

    eprintln!("Going through the BLAST hits");

    // in the original code, they save the assembly to a hash file
    // while this might be better for small genomes, this will take
    // up a huge amount of space. I'm opting for iteration but noting
    // that samtools faidx will probably be a better option

    // an iterator over the genome
    let genome_reader = fasta::Reader::from_file(genome_path.clone())?;

    for record in genome_reader.records() {
        let rec = record?;
        let query_name = rec.id();

        let outfile = blastn_dir.clone().join(format!("{}.fa", query_name));

        // create a new fasta writer with the name of the output as the
        // id name of the record. We want to save this in the blastn directory

        // the file first, in the blastn dir
        let fasta_path = blastn_dir.clone().join(format!("{}.fa", query_name));
        // now the file
        let fasta_file = fs::File::create(&fasta_path)?;
        // and the handle
        let handle = io::BufWriter::new(fasta_file);
        let mut writer = fasta::Writer::new(handle);

        // write the record to the fasta file
        writer.write_record(&rec)?;

        // now open up the blast output file
        let blast_hits = BlastRecord::from_file(temp_blast_out.clone())?;
        // filter by query name
        let mut filtered_blast_hits = blast_hits.filter_by_query_name(query_name);
        // sort blast hits on evalue, and select top 20 records
        filtered_blast_hits.sort_by_evalue();
        let top_blast_hits = filtered_blast_hits.top_n(20);
        let unique_hits = top_blast_hits.clone().filter_unique_combinations();

        for hit in unique_hits.0 {
            let query = hit.qseqid;
            let subject = hit.sseqid;

            // count minus and plus strands
            let (mut minus, mut plus) = (0, 0);

            let top_blast_hits_inner = top_blast_hits.clone();

            // Filter hits for this query-subject pair and sort by alignment positions
            let mut filtered_hits = top_blast_hits_inner.filter_by_query_subject(&query, &subject);
            filtered_hits.sort_by_alignment_positions();

            let mut start = 0;
            let mut stop = 0;

            for (index, hit) in filtered_hits.0.iter().enumerate() {
                let hit_start = hit.sstart;
                let hit_stop = hit.send;

                if index == 0 {
                    if hit_stop > hit_start {
                        plus += 1;
                    } else {
                        minus += 1;
                        start = hit_stop;
                        stop = hit_start;
                    }
                } else {
                    // all other lines
                    // should we merge with the previous line?
                    if min(hit_start, hit_stop) - max(start, stop) < maxhitdist as u64 {
                        start = min(min(start, hit_start), min(stop, hit_stop));
                        stop = max(max(start, hit_start), max(stop, hit_stop));

                        if hit_stop > hit_start {
                            plus += 1;
                        } else {
                            minus += 1;
                        }
                    } else {
                        // print previous result and save this line to compare with
                        // find the proper direction
                        let direction = if plus > minus { '+' } else { '-' };
                        // check that a majority of the hits has the same direction
                        if (max(plus, minus) as f64 / (plus + minus) as f64) < minfrac {
                            eprintln!("Ambiguous orientation of BLAST hits for {}", query);
                        }

                        // extract sequences and save somewhere
                        extract_seq(
                            genome_path.clone(),
                            &query,
                            start - flank,
                            stop + flank,
                            direction,
                            outfile.clone(),
                        )?;

                        // save the current line
                        // TODO: check this...
                        if hit_stop > hit_start {
                            plus = 1;
                            minus = 0;
                            start = hit_start;
                            stop = hit_stop;
                        } else {
                            plus = 0;
                            minus = 1;
                            start = hit_stop;
                            stop = hit_start;
                        }
                    }
                }
                // last line not printed yet
                // do it here:

                let direction = if plus > minus { '+' } else { '-' };
                // check that a majority of the hits has the same direction
                if (max(plus, minus) as f64 / (plus + minus) as f64) < minfrac {
                    eprintln!("Ambiguous orientation of BLAST hits for {}", query);
                }

                // extract sequences and save somewhere
                extract_seq(
                    genome_path.clone(),
                    &query,
                    start - flank,
                    stop + flank,
                    direction,
                    outfile.clone(),
                )?;

                // look at the current line
                if hit_stop > hit_start {
                    plus = 1;
                    minus = 0;
                    start = hit_start;
                    stop = hit_stop;
                } else {
                    plus = 0;
                    minus = 1;
                    start = hit_stop;
                    stop = hit_start;
                }
            }

            // last line not printed yet
            // do it here:
            let direction = if plus > minus { '+' } else { '-' };
            // check that a majority of the hits has the same direction
            if (max(plus, minus) as f64 / (plus + minus) as f64) < minfrac {
                eprintln!("Ambiguous orientation of BLAST hits for {}", query);
            }

            // extract sequences and save somewhere
            extract_seq(
                genome_path.clone(),
                &query,
                start - flank,
                stop + flank,
                direction,
                outfile.clone(),
            )?;
        }
    }

    Ok(())
}

// use an iterator approach to extract sequences
fn extract_seq(
    genome_path: PathBuf,
    query: &str,
    start: u64,
    stop: u64,
    direction: char,
    outfile: PathBuf,
) -> Result<()> {
    let genome_reader = fasta::Reader::from_file(genome_path)?;

    // Open the FASTA file in append mode
    let file = OpenOptions::new()
        .append(true) // Open in append mode
        .create(true) // Create if it doesn't exist
        .open(outfile)?;

    let mut writer = Writer::new(BufWriter::new(file));

    // get the sequence we need
    for record in genome_reader.records() {
        let rec = record?;
        if rec.id() == query {
            let mut substr = rec.seq()[start as usize - 1..stop as usize].to_vec();
            if direction == '-' {
                substr = bio::alphabets::dna::revcomp(substr);
            }
            writer.write(rec.id(), None, &substr)?;
        }
    }

    Ok(())
}
