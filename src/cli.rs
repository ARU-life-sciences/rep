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

fn recover_copies_command() -> clap::Command {
    use clap::{Arg, Command};
    Command::new("recover-copies")
        .about("Recover, extract and align genomic copies of candidate repeat families.")
        .arg(
            Arg::new("queries")
                .long("queries")
                .required(true)
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("genome")
                .long("genome")
                .required(true)
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("out")
                .long("out")
                .required(true)
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("flank-length")
                .long("flank-length")
                .default_value("2000")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("hits-per-family")
                .long("hits-per-family")
                .default_value("20")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("max-fragment-gap")
                .long("max-fragment-gap")
                .default_value("10000")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("orientation-threshold")
                .long("orientation-threshold")
                .default_value("0.8")
                .value_parser(value_parser!(f64)),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .default_value("1")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("search-backend")
                .long("search-backend")
                .default_value("blastn"),
        )
        .arg(
            Arg::new("blast-task")
                .long("blast-task")
                .default_value("dc-megablast")
                .value_parser(["megablast", "dc-megablast", "blastn", "blastn-short"]),
        )
        .arg(
            Arg::new("evalue")
                .long("evalue")
                .default_value("1e-10")
                .value_parser(value_parser!(f64)),
        )
        .arg(Arg::new("aligner").long("aligner").default_value("mafft"))
        .arg(
            Arg::new("diversity-filter")
                .long("diversity-filter")
                .help("Skip selected loci starting within 100 bp of an already selected locus.")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("force")
                .long("force")
                .action(ArgAction::SetTrue)
                .conflicts_with("resume"),
        )
        .arg(Arg::new("resume").long("resume").action(ArgAction::SetTrue))
}

pub fn parse_recover_copies_args() -> Result<crate::CopyRecoveryConfig> {
    // argv[0] is the binary and argv[1] is the `recover-copies` subcommand name
    recover_copies_from(std::env::args_os().skip(1))
}

fn recover_copies_from<I, T>(args: I) -> Result<crate::CopyRecoveryConfig>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let m = recover_copies_command()
        .try_get_matches_from(args)
        .unwrap_or_else(|e| e.exit());
    let path = |n: &str| m.get_one::<PathBuf>(n).cloned().expect("required by clap");
    let num = |n: &str| *m.get_one::<usize>(n).expect("defaulted by clap");
    let config = crate::CopyRecoveryConfig {
        genome: path("genome"),
        queries: path("queries"),
        output_dir: path("out"),
        flank_length: num("flank-length"),
        hits_per_family: num("hits-per-family"),
        max_fragment_gap: num("max-fragment-gap"),
        orientation_threshold: *m
            .get_one::<f64>("orientation-threshold")
            .expect("defaulted"),
        threads: num("threads"),
        search_backend: m
            .get_one::<String>("search-backend")
            .expect("defaulted")
            .parse()
            .map_err(|e: String| Error::new(ErrorKind::GenericCli(e)))?,
        blast_task: m
            .get_one::<String>("blast-task")
            .expect("defaulted")
            .clone(),
        evalue: *m.get_one::<f64>("evalue").expect("defaulted"),
        aligner: m
            .get_one::<String>("aligner")
            .expect("defaulted")
            .parse()
            .map_err(|e: String| Error::new(ErrorKind::GenericCli(e)))?,
        force: m.get_flag("force"),
        resume: m.get_flag("resume"),
        diversity_filter: m.get_flag("diversity-filter"),
        reputils: None,
    };
    for (label, p) in [("genome", &config.genome), ("queries", &config.queries)] {
        if !p.is_file() {
            return Err(Error::new(ErrorKind::GenericCli(format!(
                "--{label} is not a regular file: {}",
                p.display()
            ))));
        }
    }
    config
        .validate()
        .map_err(|e| Error::new(ErrorKind::GenericCli(e)))?;
    Ok(config)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recover_copies_parses_typed_options() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("x.fa");
        std::fs::write(&fa, ">a\nACGT\n").unwrap();
        let fa = fa.to_str().unwrap();
        let c = recover_copies_from([
            "recover-copies",
            "--queries",
            fa,
            "--genome",
            fa,
            "--out",
            "o",
            "--evalue",
            "1e-5",
            "--blast-task",
            "blastn",
            "--diversity-filter",
        ])
        .unwrap();
        assert_eq!(c.blast_task, "blastn");
        assert_eq!(c.evalue, 1e-5);
        assert!(c.diversity_filter);
        assert_eq!(c.hits_per_family, 20);
    }

    #[test]
    fn recover_copies_rejects_missing_inputs() {
        let err = recover_copies_from([
            "recover-copies",
            "--queries",
            "/nonexistent.fa",
            "--genome",
            "/nonexistent.fa",
            "--out",
            "o",
        ])
        .unwrap_err();
        assert!(err.to_string().contains("not a regular file"));
    }
}
