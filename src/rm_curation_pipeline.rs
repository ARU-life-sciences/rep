// give proper acknowledgement to the original source code and provide a link to the original source code.
use std::{borrow::Cow, fs, io, path::PathBuf};

use crate::{parse_blast::BlastRecord, CliArgs, Error, ErrorKind, Result, DATA, INTERMEDIATE};

use bio::io::fasta;

// TODO: before we start, have to make sure there
// are no special characters in the fasta headers
// `/`. Do this here.

pub fn rm_curation_pipeline(matches: CliArgs) -> Result<()> {
    let mut data_path = matches.configure.clone();
    data_path.push(DATA);

    // TODO: check this is correct. The aim here is to use the
    // fasta we copied to the data directory. So we need its name.
    let cli_matches_fasta_file = matches.fasta_file.clone();
    let genome_fasta_name = match cli_matches_fasta_file.file_name() {
        Some(n) => n.to_string_lossy(),
        None => {
            return Err(Error::new(ErrorKind::GenericCli(
                "Could not get fasta name".into(),
            )))
        }
    };

    // first of all we need to run "makeblastdb"
    // with the input fasta,
    // but we want the dbtype to be nucleotide and
    // parse the seqids
    let makeblastdb = std::process::Command::new("makeblastdb")
        .current_dir(data_path)
        .args(["-in", &genome_fasta_name])
        // this is what is written here:
        // https://github.com/ValentinaPeona/TardigraTE/blob/main/Practicals/Practical3.md
        // but we can modify the name later
        .args(["-out", &genome_fasta_name])
        .args(["-dbtype", "nucl"])
        .arg("-parse_seqids")
        .spawn()?;

    let _blastdb = makeblastdb.wait_with_output()?;

    // a random path so no errors for now
    let rmdl_library = PathBuf::new();

    run_rmdl_curation_pipeline(matches, genome_fasta_name, rmdl_library)?;

    Ok(())
}

// the actual pipeline
fn run_rmdl_curation_pipeline(
    matches: CliArgs,
    genome_fasta_name: Cow<str>,
    rmdl_library: PathBuf,
) -> Result<()> {
    // FIXME: this path is used a bunch of times, maybe we can
    // factor it out.
    // the path to the genome, and also the blast database
    let mut genome_and_blastdb_path = matches.configure.clone();
    genome_and_blastdb_path.push(DATA);
    genome_and_blastdb_path.push(genome_fasta_name.to_string());

    // # Make folders blastn and aligned
    // we want these folders to be in the intermediate directory

    let blastn_dir = matches.configure.join(INTERMEDIATE).join("blastn");
    fs::create_dir_all(blastn_dir)?;
    let aligned_dir = matches.configure.join(INTERMEDIATE).join("aligned");
    fs::create_dir_all(aligned_dir)?;

    let flank = 2000;
    let maxhitdist = 10000;
    let minfrac = 0.8;
    let digits = 5;
    let hits = 20;

    // TODO: keep in mind these will probably need to
    // be in context of the matches.configure
    let temp_blast_out = matches
        .configure
        .join(INTERMEDIATE)
        .join("tempBlastOut.txt");
    let temp_map_names = matches
        .configure
        .join(INTERMEDIATE)
        .join("tempMapNames.txt");
    let temp_mafft = matches.configure.join(INTERMEDIATE).join("tempMafft.txt");

    // 1. blast all files in the repeatmasked directory
    blast_repeatmasked(
        rmdl_library,
        genome_and_blastdb_path,
        temp_blast_out,
        temp_map_names,
        hits,
    )?;

    Ok(())
}

// 1. blast all files in the repeatmasked directory
fn blast_repeatmasked(
    rmdl_library: PathBuf,
    blastdb: PathBuf,
    temp_blast: PathBuf,
    temp_map_names: PathBuf,
    // FIXME: do I need this??
    _hits: i32,
) -> Result<()> {
    eprintln!(
        "rep :: blast_repeatmasked :: BLAST all sequences in {:?} against {:?} ",
        rmdl_library, blastdb
    );

    // remove the old temp map names file if it's present
    if temp_map_names.exists() {
        fs::remove_file(&temp_map_names)?;
    }

    // create a file with the name of the temp_map_names
    // FIXME: what's the point of this??
    // let temp_map_names_file = fs::File::create(&temp_map_names)?;

    // TODO: check this is okay not to use `current_dir()`
    // I think as all inputs and outputs have paths specified, this is
    // okay.
    let blastn_command = std::process::Command::new("blastn")
        .args(["-db", &blastdb.to_string_lossy()])
        .args(["-query", &rmdl_library.to_string_lossy()])
        // 7 includes header lines
        .args(["-outfmt", "7"])
        .arg("-evalue")
        .arg("10e-10")
        .args(["-out", &temp_blast.to_string_lossy()])
        .spawn()?;

    let _ = blastn_command.wait_with_output()?;

    Ok(())
}

// 2. Find hits from assembly (genome) and add original query
// &step2($FASTA, $blastdir, $maxhitdist, $minfrac, $ASSEMBLY);

fn find_hits_from_assembly(
    matches: CliArgs,
    rmdl_library: PathBuf,
    blastn_dir: PathBuf,
    maxhitdist: i32,
    minfrac: f64,
    genome_fasta_name: Cow<str>,
    temp_blast_out: PathBuf,
) -> Result<()> {
    // the path to the genome, and also the blast database
    let mut genome_path = matches.configure.clone();
    genome_path.push(DATA);
    genome_path.push(genome_fasta_name.to_string());

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
    let genome_reader = fasta::Reader::from_file(genome_path)?;
    // an iterator over our repeatmasked library
    let rmdl_library_reader = fasta::Reader::from_file(rmdl_library)?;

    // now we iterate over the rmdl_library_reader
    // in each iteration we:
    // 1. save each record as a separate fasta file in the blastdir
    // 2.

    for record in genome_reader.records() {
        let rec = record?;

        // create a new fasta writer with the name of the output as the
        // id name of the record. We want to save this in the blastn directory

        // the file first, in the blastn dir
        let fasta_path = blastn_dir.clone().join(format!("{}.fa", rec.id()));
        // now the file
        let fasta_file = fs::File::create(&fasta_path)?;
        // and the handle
        let handle = io::BufWriter::new(fasta_file);
        let mut writer = fasta::Writer::new(handle);

        // write the record to the fasta file
        writer.write_record(&rec)?;

        // now open up the blast output file
        let mut blast_hits = BlastRecord::from_file(temp_blast_out.clone())?;
        // sort blast hits on evalue, and select top 20 records
        blast_hits.sort_by_evalue();
        let top_blast_hits = blast_hits.top_n(20);
    }

    Ok(())
}
