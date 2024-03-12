use clap::{arg, command, value_parser};

use crate::error::{Error, ErrorKind, Result};
use std::fs;
use std::path::PathBuf;
use std::process::Command;

// a struct to contain all the cliargs
// at the moment we only want the path
// to the fasta file
#[derive(Debug)]
pub struct CliArgs {
    // path to the fasta file
    pub fasta_file: PathBuf,
    // whether to configure file system
    pub configure: PathBuf,
}

// check that we have the following
// executables:
// the perl scripts are optional at the moment
// calcDivergenceFromAlign.pl
// createRepeatLandscape.pl
// rmOut2Fasta.pl
// rmOutToGFF3.pl
fn check_executables() -> Result<()> {
    // automate the checking...
    eprintln!("Checking for required executables...");
    fn check_executables_inner(exec: String) -> Result<()> {
        match Command::new(exec.clone()).output() {
            Ok(_) => eprintln!("{} found", exec),
            Err(err) => {
                let error_kind = err.kind();

                eprintln!("{} not found", exec);
                eprintln!("Please install RepeatMasker/RepeatModeler and add it to your PATH");
                eprintln!("https://www.repeatmasker.org/");
                return Err(Error::new(ErrorKind::IO(error_kind)));
            }
        }
        Ok(())
    }

    // iterate over the executables
    // and run check_executables_inner
    for exec in vec![
        "RepeatMasker",
        "RepeatModeler",
        // "calcDivergenceFromAlign.pl",
        // "createRepeatLandscape.pl",
        // "rmOut2Fasta.pl",
        // "rmOutToGFF3.pl",
    ] {
        check_executables_inner(exec.to_string())?;
    }

    Ok(())
}

fn parse_args() -> Result<CliArgs> {
    let matches = command!()
        // not optional
        .arg(arg!(<FASTA> "Input file in fasta format").value_parser(value_parser!(PathBuf)))
        // not optional
        .arg(
            arg!(-c --configure <CONFIG_PATH> "Configure the file system")
                .required(true)
                .value_parser(value_parser!(PathBuf)),
        )
        .get_matches();

    // parse the arguments out
    let fasta = matches
        .get_one::<PathBuf>("FASTA")
        .cloned()
        .expect("defaulted by clap");

    let configure = matches
        .get_one::<PathBuf>("configure")
        .cloned()
        .expect("defaulted by clap");

    // collect the arguments
    Ok(CliArgs {
        fasta_file: fasta,
        configure,
    })
}

// the entry point for the whole program
pub fn pipeline() -> Result<()> {
    // check whether the executables are there first
    check_executables()?;

    // now parse the args
    let matches = parse_args()?;

    // set up the file system at the specified path
    set_up_filesystem(matches)?;

    Ok(())
}

// a function to set up the file system
// we want to create a set of directories
// at the matches.configure path
// three directories:
// 1. intermediate
// 2. results
// 3. pipeline_scripts
fn set_up_filesystem(matches: CliArgs) -> Result<()> {
    let mut configure = matches.configure;

    // make the configuration directory
    // and all the subdirectories
    configure.push("intermediate");
    fs::create_dir(configure.clone())?;
    configure.pop();
    configure.push("results");
    fs::create_dir(configure.clone())?;
    configure.pop();
    configure.push("pipeline_scripts");
    fs::create_dir(configure.clone())?;

    Ok(())
}
