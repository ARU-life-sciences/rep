use clap::{arg, command, value_parser, ArgAction};

use crate::{error::Result, Error, ErrorKind};
use std::path::{Path, PathBuf};

// a struct to contain all the CliArgs
// at the moment we only want the path
// to the fasta file
#[derive(Debug, Clone)]
pub struct CliArgs {
    // path to the primary FASTA file
    pub fasta_file: PathBuf,
    // optional separate genome searched during consensus curation
    pub curation_target: Option<PathBuf>,
    // whether to configure file system
    pub configure: Option<PathBuf>,
    // the name of the database for BuildDatabase
    pub database: Option<String>,
    // repeat modeler threads
    pub rmo_threads: u8,
    // repeat masker threads
    pub rma_threads: u8,
    // run repeat masker only
    pub rma_only: bool,
    // run consensus curation from existing RepeatModeler output
    pub curation_only: bool,
    // verbose flag
    pub verbose: bool,
}

impl CliArgs {
    pub fn configure_path(&self) -> Result<&Path> {
        self.configure.as_deref().ok_or_else(|| {
            Error::new(ErrorKind::GenericCli(
                "A configuration directory is required".into(),
            ))
        })
    }

    pub fn database_name(&self) -> Result<&str> {
        self.database.as_deref().ok_or_else(|| {
            Error::new(ErrorKind::GenericCli(
                "A database name is required to run RepeatModeler".into(),
            ))
        })
    }
}

pub fn parse_args() -> Result<CliArgs> {
    let matches = command!()
        .next_line_help(true)
        // not optional
        .arg(
            arg!(<FASTA> "Input file in FASTA format.")
                .value_parser(value_parser!(PathBuf)
                )
        )
        .arg(
            arg!(--curation_target <CURATION_FASTA> "Optional separate FASTA searched for consensus copies.")
                .value_parser(value_parser!(PathBuf)),
        )
        // not optional
        .arg(
            arg!(-c --configure <CONFIG_PATH> "Configure the file system - and create the required directories.")
                .required(true)
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            arg!(-d --database <DATABASE_NAME> "Name of the database, when building using `BuildDatabase`.")
                // not required unless you only want to run RepeatMasker
                // or I think, when you want run the curation pipeline
                .required_unless_present_any(["rma_only", "curation_only"])
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
            arg!(--rma_only "Run RepeatMasker only. Skip RepeatModeler and curation.")
                .conflicts_with("curation_only")
                .action(ArgAction::SetTrue)
        )
        .arg(
            arg!(--curation_only "Run curation from existing RepeatModeler output.")
                .action(ArgAction::SetTrue)
        )
        .arg(arg!(--verbose "Print extra debug information").action(ArgAction::SetTrue))
        .get_matches();

    // parse the arguments out
    let fasta = matches
        .get_one::<PathBuf>("FASTA")
        .cloned()
        .expect("errored by clap");

    if !fasta.is_file() {
        return Err(Error::new(ErrorKind::GenericCli(format!(
            "FASTA path is not a regular file: {}",
            fasta.display()
        ))));
    }

    let curation_target = matches.get_one::<PathBuf>("curation_target").cloned();
    if let Some(path) = &curation_target {
        if !path.is_file() {
            return Err(Error::new(ErrorKind::GenericCli(format!(
                "Curation target is not a regular file: {}",
                path.display()
            ))));
        }
    }

    let configure = matches.get_one::<PathBuf>("configure").cloned();

    let database = matches.get_one::<String>("database").cloned();

    let rmo_threads = matches
        .get_one::<u8>("rmo_threads")
        .cloned()
        .expect("errored by clap");

    let rma_threads = matches
        .get_one::<u8>("rma_threads")
        .cloned()
        .expect("errored by clap");

    let rma_only = matches.get_flag("rma_only");

    let curation_only = matches.get_flag("curation_only");

    let verbose = matches.get_flag("verbose");

    // collect the arguments
    Ok(CliArgs {
        fasta_file: fasta,
        curation_target,
        configure,
        database,
        rmo_threads,
        rma_threads,
        rma_only,
        curation_only,
        verbose,
    })
}

pub fn parse_recover_copies_args() -> Result<crate::CopyRecoveryConfig> {
    use clap::{Arg, Command};
    let m = Command::new("recover-copies")
        .arg(Arg::new("queries").long("queries").required(true))
        .arg(Arg::new("genome").long("genome").required(true))
        .arg(Arg::new("out").long("out").required(true))
        .arg(
            Arg::new("flank-length")
                .long("flank-length")
                .default_value("2000"),
        )
        .arg(
            Arg::new("hits-per-family")
                .long("hits-per-family")
                .default_value("20"),
        )
        .arg(
            Arg::new("max-fragment-gap")
                .long("max-fragment-gap")
                .default_value("10000"),
        )
        .arg(
            Arg::new("orientation-threshold")
                .long("orientation-threshold")
                .default_value("0.8"),
        )
        .arg(Arg::new("threads").long("threads").default_value("1"))
        .arg(
            Arg::new("search-backend")
                .long("search-backend")
                .default_value("blastn"),
        )
        .arg(Arg::new("aligner").long("aligner").default_value("mafft"))
        .arg(Arg::new("force").long("force").action(ArgAction::SetTrue))
        .arg(Arg::new("resume").long("resume").action(ArgAction::SetTrue))
        .get_matches_from(std::env::args().skip(1));
    let get = |n: &str| m.get_one::<String>(n).unwrap();
    Ok(crate::CopyRecoveryConfig {
        genome: get("genome").into(),
        queries: get("queries").into(),
        output_dir: get("out").into(),
        flank_length: get("flank-length").parse().unwrap(),
        hits_per_family: get("hits-per-family").parse().unwrap(),
        max_fragment_gap: get("max-fragment-gap").parse().unwrap(),
        orientation_threshold: get("orientation-threshold").parse().unwrap(),
        threads: get("threads").parse().unwrap(),
        search_backend: get("search-backend")
            .parse()
            .map_err(|e: String| Error::new(ErrorKind::GenericCli(e)))?,
        aligner: get("aligner")
            .parse()
            .map_err(|e: String| Error::new(ErrorKind::GenericCli(e)))?,
        force: m.get_flag("force"),
        resume: m.get_flag("resume"),
        diversity_filter: false,
        reputils: None,
    })
}
