// the error module
pub mod error;

// where all cli related stuff goes
pub mod cli;

pub use cli::{check_executables, parse_args, CliArgs};
pub use error::Result;
use std::fs;

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
    fs::create_dir_all(configure.clone())?;
    configure.pop();
    configure.push("results");
    fs::create_dir_all(configure.clone())?;
    configure.pop();
    configure.push("pipeline_scripts");
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

    Ok(())
}
