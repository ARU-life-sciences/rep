// the error module
pub mod error;

// where all cli related stuff goes
pub mod cli;

// where we will run RepeatModeler
pub mod repeatmodeler;

pub use cli::{parse_args, CliArgs};
pub use error::{Error, ErrorKind, Result};
pub use repeatmodeler::run_repeatmodeler;
use std::{
    fs::{self, File},
    process::{Command, Stdio},
};

// the output paths here
const INTERMEDIATE: &str = "intermediate";
const RESULTS: &str = "results";
const PIPELINE_SCRIPTS: &str = "pipeline_scripts";
const DATA: &str = "data";

// the entry point for the whole program
pub fn pipeline() -> Result<()> {
    // check whether the executables are there first
    check_executables()?;

    // now parse the args
    let matches = parse_args()?;

    // set up the file system at the specified path
    set_up_filesystem(matches.clone())?;

    // and now we need to actually run the analyses.
    run_repeatmodeler(matches)?;

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

// a function to set up the file system
// we want to create a set of directories
// at the matches.configure path
// three directories:
// 1. intermediate
// 2. results
// 3. pipeline_scripts
// 4. data
fn set_up_filesystem(matches: CliArgs) -> Result<()> {
    let mut configure = matches.configure;

    // make the configuration directory
    // and all the subdirectories
    configure.push(INTERMEDIATE);
    fs::create_dir_all(configure.clone())?;
    configure.pop();
    configure.push(RESULTS);
    fs::create_dir_all(configure.clone())?;
    configure.pop();
    configure.push(PIPELINE_SCRIPTS);
    fs::create_dir_all(configure.clone())?;

    // also push the code from the `perl` folder
    // into the pipeline_scripts folder
    let rmdl_curation_pipeline = include_str!("perl/RMDL_curation_pipeline.pl");
    let rename_rmdl_consensi = include_str!("perl/renameRMDLconsensi.pl");
    let shorten_scaffold_names = include_str!("perl/shortenScaffoldnames.pl");

    // now write these to file
    for (code, path) in [
        (rmdl_curation_pipeline, "RMDL_curation_pipeline.pl"),
        (rename_rmdl_consensi, "renameRMDLconsensi.pl"),
        (shorten_scaffold_names, "shortenScaffoldnames.pl"),
    ]
    .iter()
    {
        configure.push(path);
        fs::write(configure.clone(), code)?;
        configure.pop();
    }

    // now deal with the data
    configure.pop();
    configure.push(DATA);
    fs::create_dir_all(configure.clone())?;

    // check the ending of the file.
    match matches.fasta_file.to_string_lossy().ends_with("gz") {
        true => {
            // the final component of the fasta file name here
            let base_fasta_name = matches.fasta_file.file_name().unwrap();
            let fasta_name = base_fasta_name.to_string_lossy();
            // FIXME: the following line will panic if the file is mal-formatted
            let new_fasta_name = fasta_name.strip_suffix(".gz").unwrap();

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

            std::io::copy(&mut gunzip_process.stdout.unwrap(), &mut f)?;
        }
        false => {
            // else use cp
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

    eprintln!(
        "Successfully copied {}",
        matches.fasta_file.to_string_lossy()
    );

    Ok(())
}
