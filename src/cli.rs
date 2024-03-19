use clap::{arg, command, value_parser};

use crate::error::Result;
use std::path::PathBuf;

// a struct to contain all the cliargs
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
    // repeat modeller threads
    pub rm_threads: u8,
}

pub fn parse_args() -> Result<CliArgs> {
    let matches = command!()
        // not optional
        .arg(arg!(<FASTA> "Input file in fasta format").value_parser(value_parser!(PathBuf)))
        // not optional
        .arg(
            arg!(-c --configure <CONFIG_PATH> "Configure the file system - and create the required directories")
                .required(true)
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            arg!(-d --database <DATABASE_NAME> "Name of the database, when building using `BuildDatabase`.")
                .required(true)
                .value_parser(value_parser!(String)),
        )
        .arg(
            arg!(--rm_threads <RM_THREADS> "Number of threads to use for RepeatModeller")
                .default_value("8")
                .value_parser(value_parser!(u8)),

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

    let rm_threads = matches
        .get_one::<u8>("rm_threads")
        .cloned()
        .expect("errored by clap");

    // collect the arguments
    Ok(CliArgs {
        fasta_file: fasta,
        configure,
        database,
        rm_threads,
    })
}
