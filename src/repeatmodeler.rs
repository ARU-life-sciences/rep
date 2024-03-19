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

    // TODO: eventually have the stderr/stdout of commands
    // as optional, maybe behind a verbose flag, and saved to a file

    // build the database here
    let build_database = Command::new("BuildDatabase")
        .current_dir(data_path.clone())
        .arg("-name")
        .arg(matches.database.clone())
        .arg("-dir")
        .arg(".")
        .spawn()?;

    let _build_database_out = build_database.wait_with_output()?;

    let run_repeat_modeler = Command::new("RepeatModeler")
        .current_dir(data_path)
        .arg("-database")
        .arg(matches.database)
        .arg("-threads")
        .arg(matches.rm_threads.to_string())
        .spawn()?;

    let _run_repeat_modeler_out = run_repeat_modeler.wait_with_output()?;

    // everything seems to be in the 'data' directory now
    // so we will have to move things around afterwards.
    Ok(())
}
