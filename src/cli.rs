use clap::{arg, command, value_parser, ArgAction};

use crate::error::Result;
use std::path::PathBuf;

// a struct to contain all the CliArgs
// at the moment we only want the path
// to the fasta file
#[derive(Debug, Clone)]
pub struct CliArgs {
    // path to the fasta file
    pub fasta_file: PathBuf,
    // whether to configure file system
    pub configure: PathBuf,
    // the name of the database for BuildDatabase
    pub database: String,
    // repeat modeler threads
    pub rmo_threads: u8,
    // repeat masker threads
    pub rma_threads: u8,
    // run repeat masker only
    pub rma_only: bool,
}

pub fn parse_args() -> Result<CliArgs> {
    let matches = command!()
        .next_line_help(true)
        // not optional
        .arg(arg!(<FASTA> "Input file in fasta format").value_parser(value_parser!(PathBuf)))
        // not optional
        .arg(
            arg!(-c --configure <CONFIG_PATH> "Configure the file system - and create the required directories.")
                // not required unless you only want to run RepeatMasker
                .required_unless_present("rma_only")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            arg!(-d --database <DATABASE_NAME> "Name of the database, when building using `BuildDatabase`.")
                // not required unless you only want to run RepeatMasker
                .required_unless_present("rma_only")
                .value_parser(value_parser!(String)),
        )
        .arg(
            arg!(--rmo_threads <RMO_THREADS> "Number of threads to use for RepeatModeler.")
                .default_value("8")
                .value_parser(value_parser!(u8)),

        )
        .arg(
            arg!(--rma_threads <RMA_THREADS> "Number of threads to use for RepeatMasker.")
                .default_value("8")
                .value_parser(value_parser!(u8)),

        )
        .arg(
            arg!(--rma_only "Run RepeatMasker only. Skip RepeatModeler; currently for development.")
                .action(ArgAction::SetTrue)
        )
        .get_matches();

    // parse the arguments out
    let fasta = matches
        .get_one::<PathBuf>("FASTA")
        .cloned()
        .expect("errored by clap");

    let configure = matches
        .get_one::<PathBuf>("configure")
        .cloned()
        .expect("errored by clap");

    let database = matches
        .get_one::<String>("database")
        .cloned()
        .expect("errored by clap");

    let rmo_threads = matches
        .get_one::<u8>("rmo_threads")
        .cloned()
        .expect("errored by clap");

    let rma_threads = matches
        .get_one::<u8>("rma_threads")
        .cloned()
        .expect("errored by clap");

    let rma_only = matches.get_flag("rma_only");

    // collect the arguments
    Ok(CliArgs {
        fasta_file: fasta,
        configure,
        database,
        rmo_threads,
        rma_threads,
        rma_only,
    })
}
