// There are two steps here
// 1. create the database for the genome using BuildDatabase
// 2. run RepeatModeler on the database

use crate::{CliArgs, Result, DATA};
use std::process::Command;

pub fn run_repeatmodeler(matches: CliArgs) -> Result<()> {
    // we have all of our directories set up.
    // we need to specify the directory with the data in it
    let mut data_path = matches.configure.clone();
    data_path.push(DATA);

    // build the database here
    let build_database = Command::new("BuildDatabase")
        .arg("-name")
        .arg(matches.database)
        .arg("-dir")
        .arg(data_path)
        .spawn()?;

    let build_database_out = build_database.wait_with_output()?;
    // TODO: eventually have the output as optional, maybe behind
    // a verbose flag
    eprintln!(
        "STDOUT: {}",
        String::from_utf8(build_database_out.stdout).unwrap()
    );
    eprintln!(
        "STDERR: {}",
        String::from_utf8(build_database_out.stderr).unwrap()
    );

    Ok(())
}
