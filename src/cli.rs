use clap::{arg, command, value_parser, Command};

use crate::error::Result;
use std::path::PathBuf;

// a struct to contain all the cliargs
// at the moment we only want the path
// to the fasta file
#[derive(Debug, Clone)]
pub struct RepeatModelerCliArgs {
    // path to the RepeatModeler executable
    pub repeatmodeler: PathBuf,
    // path to the fasta file
    pub fasta_file: PathBuf,
    // whether to configure file system
    pub configure: PathBuf,
    // the name of the database for BuildDatabase
    pub database: String,
    // repeat modeller threads
    pub rm_threads: u8,
}

#[derive(Debug, Clone)]
pub struct RepeatMaskerCliArgs {
    // the repeatmasker executable
    pub repeatmasker: PathBuf,
    // the classified repeat library
    pub consensi_classified: PathBuf,
    // the configured directory
    pub configure: PathBuf,
    // the directory to output the results
    pub dir: PathBuf,
    // the genome to mask/identify repeats in
    pub genome: PathBuf,
    // number of threads
    pub rm_threads: u8,
}

// for the subcommands
#[derive(Debug, Clone)]
pub enum CliArgs {
    RepeatModeler(RepeatModelerCliArgs),
    RepeatMasker(RepeatMaskerCliArgs),
}

pub fn parse_args() -> Result<CliArgs> {
    let matches = command!()
        .propagate_version(true)
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(
            Command::new("repeatmodeler")
                .about("Build a RepeatModeler database and run RepeatModeler")
                .arg(arg!(-r --repeatmodeler_path [RMPATH] "The path to the RepeatModeler executable")
                        .default_value("RepeatModeler")
                        .value_parser(value_parser!(PathBuf))
                )
                // not optional
                .arg(arg!(<FASTA> "Input file in fasta format").value_parser(value_parser!(PathBuf)))
                // not optional
                .arg(
                    arg!(-c --configure <CONFIG_PATH> "Configure the file system - and create the required directories")
                        .required(true)
                        .value_parser(value_parser!(PathBuf)),
                )
                .arg(
                    arg!(-d --database <DATABASE_NAME> "Name of the database, when building using `BuildDatabase`")
                        .required(true)
                        .value_parser(value_parser!(String)),
                )
                .arg(
                    arg!(--rm_threads <RM_THREADS> "Number of threads to use for RepeatModeler")
                        .default_value("8")
                        .value_parser(value_parser!(u8)),
                )

        )
        .subcommand(Command::new("repeatmasker").about("Run RepeatMasker")
            .arg(arg!(-r --repeatmasker_path [RMPATH] "The path to the RepeatMasker executable")
                .default_value("RepeatMasker")
                .value_parser(value_parser!(PathBuf))
            )
            .arg(arg!(-c --consensi_classified <CONSENSI_CLASSIFIED> "The classified repeat library")
                .required(true)
                .value_parser(value_parser!(PathBuf))
            )
            .arg(arg!(-d --dir <DIR> "The directory to output the results")
                .required(true)
                .value_parser(value_parser!(PathBuf))
            )
            .arg(arg!(<GENOME> "The genome to mask/identify repeats in").value_parser(value_parser!(PathBuf)))
            .arg(arg!(--rm_threads <RM_THREADS> "Number of threads to use for RepeatMasker")
                .default_value("8")
                .value_parser(value_parser!(u8))
            )
        )
        .get_matches();

    match matches.subcommand() {
        Some(("repeatmodeler", repeatmodeler_matches)) => {
            // parse the arguments out
            let repeatmodeler = repeatmodeler_matches
                .get_one::<PathBuf>("repeatmodeler_path")
                .cloned()
                .expect("errored by clap");

            let fasta = repeatmodeler_matches
                .get_one::<PathBuf>("FASTA")
                .cloned()
                .expect("errored by clap");

            let configure = repeatmodeler_matches
                .get_one::<PathBuf>("configure")
                .cloned()
                .expect("errored by clap");

            let database = repeatmodeler_matches
                .get_one::<String>("database")
                .cloned()
                .expect("errored by clap");

            let rm_threads = repeatmodeler_matches
                .get_one::<u8>("rm_threads")
                .cloned()
                .expect("errored by clap");

            // collect the arguments
            Ok(CliArgs::RepeatModeler(RepeatModelerCliArgs {
                repeatmodeler,
                fasta_file: fasta,
                configure,
                database,
                rm_threads,
            }))
        }
        Some(("repeatmasker", repeatmodeler_matches)) => {
            let repeatmasker = repeatmodeler_matches
                .get_one::<PathBuf>("repeatmasker_path")
                .cloned()
                .expect("errored by clap");

            let consensi_classified = repeatmodeler_matches
                .get_one::<PathBuf>("consensi_classified")
                .cloned()
                .expect("errored by clap");

            let configure = repeatmodeler_matches
                .get_one::<PathBuf>("configure")
                .cloned()
                .expect("errored by clap");

            let dir = repeatmodeler_matches
                .get_one::<PathBuf>("dir")
                .cloned()
                .expect("errored by clap");

            let genome = repeatmodeler_matches
                .get_one::<PathBuf>("GENOME")
                .cloned()
                .expect("errored by clap");

            let rm_threads = repeatmodeler_matches
                .get_one::<u8>("rm_threads")
                .cloned()
                .expect("errored by clap");

            Ok(CliArgs::RepeatMasker(RepeatMaskerCliArgs {
                repeatmasker,
                consensi_classified,
                configure,
                dir,
                genome,
                rm_threads,
            }))
        }
        _ => unreachable!(),
    }
}
