// Public modules used across the CLI tool
pub mod cli; // Command-line argument parsing
pub mod command_runner;
pub mod error; // Error types and handling
pub mod parse_blast; // BLAST outfmt 7 parser
pub mod repeatmasker; // RepeatMasker wrapper
pub mod repeatmodeler; // RepeatModeler wrapper

// Re-export key types and functions
pub use cli::{parse_args, CliArgs};
pub use command_runner::{CommandRunner, RealCommandRunner};
pub use error::{Error, ErrorKind, Result};
pub use repeatmasker::run_repeatmasker;
pub use repeatmodeler::run_repeatmodeler;

use std::{
    fs::{self, File},
    path::PathBuf,
    process::{Command, Stdio},
};

// subdirectories used in the pipeline
const INTERMEDIATE: &str = "intermediate";
const RESULTS: &str = "results";
const PIPELINE_SCRIPTS: &str = "pipeline_scripts";
const DATA: &str = "data";

// Utility to create a named subdirectory within a base path
fn make_subdir(base: &PathBuf, name: &str) -> Result<()> {
    let mut p = base.clone();
    p.push(name);
    fs::create_dir_all(&p)?;
    Ok(())
}

// Main entry point for running the full pipeline
pub fn pipeline() -> Result<()> {
    // now parse the args
    let matches = parse_args()?;

    // check whether the executables are there first
    check_executables()?;

    // set up the file system at the specified path
    set_up_filesystem(matches.clone())?;

    if matches.rma_only {
        // if we are just running repeatmodeler
        // then run it and exit
        eprintln!("Running RepeatMasker only...");
        let runner = RealCommandRunner;
        run_repeatmasker(matches.clone(), &runner)?;
        return Ok(());
    }

    // and now we need to actually run the analyses.
    eprintln!("Running RepeatModeler...");
    let runner = RealCommandRunner;
    run_repeatmodeler(matches.clone(), &runner)?;

    // and also run repeatmasker
    eprintln!("Running RepeatMasker...");
    run_repeatmasker(matches.clone(), &runner)?;

    Ok(())
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

                // TODO: move this printing to the error module
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
    for exec in [
        "/software/team301/repeat-annotation/RepeatMasker/RepeatMasker",
        "/software/team301/repeat-annotation/RepeatModeler-2.0.5/RepeatModeler",
    ] {
        check_executables_inner(exec.to_string())?;
    }

    Ok(())
}

// a function to set up the file system
// we want to create a set of directories
// at the matches.configure path
// three directories:
// 1. intermediate
// 2. results
// 3. pipeline_scripts
// 4. data
//   - RepeatModeler data
//   - RepeatMasker data
fn set_up_filesystem(matches: CliArgs) -> Result<()> {
    let mut configure = matches.configure.clone().unwrap();

    // make the configuration directory
    // and all the subdirectories
    make_subdir(&configure, INTERMEDIATE)?;
    make_subdir(&configure, RESULTS)?;
    make_subdir(&configure, PIPELINE_SCRIPTS)?;
    make_subdir(&configure, DATA)?;

    configure.push(DATA);

    // check the ending of the file.
    match matches.fasta_file.to_string_lossy().ends_with("gz") {
        true => {
            // the final component of the fasta file name here
            let base_fasta_name = matches.fasta_file.file_name().unwrap();
            let fasta_name = base_fasta_name.to_string_lossy();
            // FIXME: the following line will panic if the file is mal-formatted
            let new_fasta_name = fasta_name.strip_suffix(".gz").ok_or_else(|| {
                Error::new(ErrorKind::GenericCli(
                    "FASTA filename does not end with .gz".into(),
                ))
            })?;

            configure.push(new_fasta_name);
            // add the new command here
            let gunzip_process = Command::new("gunzip")
                .arg("-c")
                .arg(matches.fasta_file.clone())
                .stdout(Stdio::piped())
                .spawn()?;

            // little bit of wizardry here
            // see https://stackoverflow.com/questions/43949612/redirect-output-of-child-process-spawned-from-rust
            // as > cannot be used in Command.
            let mut f = File::create(&configure)?;

            // FIXME: remove this unwrap
            std::io::copy(&mut gunzip_process.stdout.unwrap(), &mut f)?;
        }
        false => {
            // else use cp
            // FIXME: remove this unwrap
            let base_fasta_name = matches.fasta_file.file_name().unwrap();
            configure.push(base_fasta_name);
            Command::new("cp")
                .arg("-r")
                .arg(matches.fasta_file.clone())
                .arg(&configure)
                .spawn()?
                .wait_with_output()?;
        }
    }

    // make separate subdir for RepeatMasker and RepeatModeler
    // within the data directory.
    eprintln!("Making repeatmasker and repeatmodeler dirs...");
    // do it again to make sure
    let mut configure = matches.configure.clone().unwrap();
    configure.push(DATA);
    configure.push("repeatmasker");
    fs::create_dir_all(configure.clone())?;
    configure.pop();
    configure.push("repeatmodeler");
    fs::create_dir_all(configure.clone())?;
    configure.pop();
    eprintln!("Made repeatmasker and repeatmodeler dirs...");

    eprintln!(
        "Successfully copied {}",
        matches.fasta_file.to_string_lossy()
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_directory_structure_created() {
        let dir = tempdir().unwrap();
        let config_path = dir.path().to_path_buf();
        let fasta_path = dir.path().join("genome.fa");
        std::fs::write(&fasta_path, ">seq\nACGT").unwrap();

        let args = CliArgs {
            fasta_file: fasta_path,
            configure: Some(config_path.clone()),
            database: Some("genomedb".into()),
            rmo_threads: 1,
            rma_threads: 1,
            rma_only: false,
            verbose: false,
        };

        set_up_filesystem(args).unwrap();

        assert!(config_path.join("data").exists());
        assert!(config_path.join("intermediate").exists());
        assert!(config_path.join("results").exists());
        assert!(config_path.join("pipeline_scripts").exists());
        assert!(config_path.join("data").join("repeatmasker").exists());
        assert!(config_path.join("data").join("repeatmodeler").exists());
    }
}
