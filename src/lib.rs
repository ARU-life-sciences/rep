// Public modules used across the CLI tool
pub mod cli; // Command-line argument parsing
pub mod command_runner;
pub mod copy_recovery;
pub mod curation;
pub mod error; // Error types and handling
pub mod parse_blast; // BLAST outfmt 7 parser
pub mod repeatmasker; // RepeatMasker wrapper
pub mod repeatmodeler; // RepeatModeler wrapper

// Re-export key types and functions
pub use cli::{parse_args, CliArgs};
pub use command_runner::{CommandRunner, RealCommandRunner};
pub use copy_recovery::{
    run_copy_recovery, Aligner, CopyRecoveryConfig, CopyRecoveryResult, FamilyRecoveryResult,
    SearchBackend,
};
pub use curation::{run_curation, run_curation_with_config, CurationConfig};
pub use error::{Error, ErrorKind, Result};
pub use repeatmasker::run_repeatmasker;
pub use repeatmodeler::run_repeatmodeler;

use std::{
    env,
    ffi::{OsStr, OsString},
    fs::{self, File},
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

// subdirectories used in the pipeline
const INTERMEDIATE: &str = "intermediate";
const RESULTS: &str = "results";
const PIPELINE_SCRIPTS: &str = "pipeline_scripts";
const DATA: &str = "data";
const REPEATMASKER_ENV: &str = "REP_REPEATMASKER";
const REPEATMODELER_ENV: &str = "REP_REPEATMODELER";
const BUILD_DATABASE_ENV: &str = "REP_BUILD_DATABASE";
const MAKEBLASTDB_ENV: &str = "REP_MAKEBLASTDB";
const BLASTN_ENV: &str = "REP_BLASTN";
const MAFFT_ENV: &str = "REP_MAFFT";

// Utility to create a named subdirectory within a base path
fn make_subdir(base: &Path, name: &str) -> Result<()> {
    let mut p = base.to_path_buf();
    p.push(name);
    fs::create_dir_all(&p)?;
    Ok(())
}

// Main entry point for running the full pipeline
pub fn pipeline() -> Result<()> {
    // now parse the args
    let mut matches = parse_args()?;

    // check whether the executables are there first
    check_executables(&matches)?;

    if matches.curation_only {
        eprintln!("Running consensus curation only...");
        let runner = RealCommandRunner;
        curation::run_curation_via_copy_recovery(matches, &runner)?;
        return Ok(());
    }

    // set up the file system at the specified path
    matches.fasta_file = set_up_filesystem(matches.clone())?;

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

    eprintln!("Running consensus curation...");
    curation::run_curation_via_copy_recovery(matches.clone(), &runner)?;

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
fn check_executables(matches: &CliArgs) -> Result<()> {
    eprintln!("Checking for required executables...");

    if matches.curation_only {
        for (env_var, executable) in [
            (MAKEBLASTDB_ENV, "makeblastdb"),
            (BLASTN_ENV, "blastn"),
            (MAFFT_ENV, "mafft"),
        ] {
            let path = find_executable(env_var, executable)?;
            eprintln!("{} found at {}", executable, path.display());
        }
        return Ok(());
    }

    let mut required = vec![(REPEATMASKER_ENV, "RepeatMasker")];
    if !matches.rma_only {
        required.push((REPEATMODELER_ENV, "RepeatModeler"));
        required.push((BUILD_DATABASE_ENV, "BuildDatabase"));
        required.push((MAKEBLASTDB_ENV, "makeblastdb"));
        required.push((BLASTN_ENV, "blastn"));
        required.push((MAFFT_ENV, "mafft"));
    }

    for (env_var, executable) in required {
        let path = find_executable(env_var, executable)?;
        eprintln!("{} found at {}", executable, path.display());
    }
    Ok(())
}

fn configured_executable(env_var: &str, default: &str) -> OsString {
    env::var_os(env_var).unwrap_or_else(|| OsString::from(default))
}

fn find_executable(env_var: &str, default: &str) -> Result<PathBuf> {
    let configured = configured_executable(env_var, default);
    let configured_path = Path::new(&configured);

    if configured_path.components().count() > 1 {
        if is_executable(configured_path) {
            return Ok(configured_path.to_path_buf());
        }
        return Err(missing_executable(env_var, &configured));
    }

    if let Some(path) = env::var_os("PATH") {
        for directory in env::split_paths(&path) {
            let candidate = directory.join(&configured);
            if is_executable(&candidate) {
                return Ok(candidate);
            }
        }
    }

    Err(missing_executable(env_var, &configured))
}

#[cfg(unix)]
fn is_executable(path: &Path) -> bool {
    use std::os::unix::fs::PermissionsExt;
    path.metadata()
        .map(|metadata| metadata.is_file() && metadata.permissions().mode() & 0o111 != 0)
        .unwrap_or(false)
}

#[cfg(not(unix))]
fn is_executable(path: &Path) -> bool {
    path.is_file()
}

fn missing_executable(env_var: &str, executable: &OsStr) -> Error {
    Error::new(ErrorKind::GenericCli(format!(
        "Required executable '{}' was not found or is not executable; add it to PATH or set {env_var}",
        executable.to_string_lossy()
    )))
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
fn set_up_filesystem(matches: CliArgs) -> Result<PathBuf> {
    let mut configure = matches.configure_path()?.to_path_buf();

    // make the configuration directory
    // and all the subdirectories
    make_subdir(&configure, INTERMEDIATE)?;
    make_subdir(&configure, RESULTS)?;
    make_subdir(&configure, PIPELINE_SCRIPTS)?;
    make_subdir(&configure, DATA)?;

    configure.push(DATA);

    // check the ending of the file.
    match matches.fasta_file.extension().and_then(|ext| ext.to_str()) == Some("gz") {
        true => {
            // the final component of the fasta file name here
            let base_fasta_name = matches.fasta_file.file_name().ok_or_else(|| {
                Error::new(ErrorKind::GenericCli("FASTA path has no filename".into()))
            })?;
            let fasta_name = base_fasta_name.to_string_lossy();
            // FIXME: the following line will panic if the file is mal-formatted
            let new_fasta_name = fasta_name.strip_suffix(".gz").ok_or_else(|| {
                Error::new(ErrorKind::GenericCli(
                    "FASTA filename does not end with .gz".into(),
                ))
            })?;

            configure.push(new_fasta_name);
            // add the new command here
            let mut gunzip_process = Command::new("gunzip")
                .arg("-c")
                .arg(matches.fasta_file.clone())
                .stdout(Stdio::piped())
                .spawn()?;

            // little bit of wizardry here
            // see https://stackoverflow.com/questions/43949612/redirect-output-of-child-process-spawned-from-rust
            // as > cannot be used in Command.
            let mut f = File::create(&configure)?;

            let mut stdout = gunzip_process.stdout.take().ok_or_else(|| {
                Error::new(ErrorKind::GenericCli(
                    "Could not capture gunzip output".into(),
                ))
            })?;
            std::io::copy(&mut stdout, &mut f)?;
            let status = gunzip_process.wait()?;
            if !status.success() {
                return Err(Error::new(ErrorKind::GenericCli(format!(
                    "gunzip failed with status {status}"
                ))));
            }
        }
        false => {
            let base_fasta_name = matches.fasta_file.file_name().ok_or_else(|| {
                Error::new(ErrorKind::GenericCli("FASTA path has no filename".into()))
            })?;
            configure.push(base_fasta_name);
            fs::copy(&matches.fasta_file, &configure)?;
        }
    }

    let staged_fasta = configure.clone();

    // make separate subdir for RepeatMasker and RepeatModeler
    // within the data directory.
    eprintln!("Making repeatmasker and repeatmodeler dirs...");
    // do it again to make sure
    let mut configure = matches.configure_path()?.to_path_buf();
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

    Ok(staged_fasta)
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
            curation_target: None,
            configure: Some(config_path.clone()),
            database: Some("genomedb".into()),
            rmo_threads: 1,
            rma_threads: 1,
            rma_only: false,
            curation_only: false,
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
