// There are two steps here
// 1. create the database for the genome using BuildDatabase
// 2. run RepeatModeler on the database

use crate::{
    command_runner::run_checked, configured_executable, CliArgs, CommandRunner, Result,
    BUILD_DATABASE_ENV, DATA, REPEATMODELER_ENV,
};
use std::process::Command;

pub fn run_repeatmodeler(matches: CliArgs, runner: &dyn CommandRunner) -> Result<()> {
    // if we are running repeatmasker only
    if matches.rma_only {
        eprintln!("Only running RepeatMasker, skipping RepeatModeler");
        return Ok(());
    }

    // we have all of our directories set up.
    // we need to specify the directory with the data in it
    let mut data_path = matches.configure_path()?.to_path_buf();
    let database = matches.database_name()?.to_owned();
    data_path.push(DATA);
    // and go into the repeatmodeler dir
    data_path.push("repeatmodeler");

    // TODO: eventually have the stderr/stdout of commands
    // as optional, maybe behind a verbose flag, and saved to a file

    // build the database here
    let mut build_database =
        Command::new(configured_executable(BUILD_DATABASE_ENV, "BuildDatabase"));

    eprintln!("Data path: {:?}", data_path);

    build_database
        .current_dir(data_path.clone())
        .arg("-name")
        .arg(&database)
        .arg("-dir")
        .arg(".")
        .arg(matches.fasta_file);

    run_checked(runner, &mut build_database, "BuildDatabase")?;

    let mut run_repeat_modeler =
        Command::new(configured_executable(REPEATMODELER_ENV, "RepeatModeler"));

    run_repeat_modeler
        .current_dir(data_path)
        .arg("-database")
        .arg(&database)
        .arg("-threads")
        .arg(matches.rmo_threads.to_string());

    run_checked(runner, &mut run_repeat_modeler, "RepeatModeler")?;

    // everything seems to be in the 'data' directory now
    // so we will have to move things around afterwards.
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::set_up_filesystem;

    use super::*;
    use std::os::unix::process::ExitStatusExt;
    use std::process::{ExitStatus, Output};
    use tempfile::tempdir;

    struct MockRunner;

    impl CommandRunner for MockRunner {
        fn run(&self, cmd: &mut Command) -> Result<Output> {
            let args = cmd
                .get_args()
                .map(|s| s.to_string_lossy().to_string())
                .collect::<Vec<_>>();
            let program = cmd.get_program().to_string_lossy();
            eprintln!("[MOCK] {} {:?}", program, args);
            Ok(Output {
                status: ExitStatus::from_raw(0),
                stdout: "mock stdout".as_bytes().to_vec(),
                stderr: "mock stderr".as_bytes().to_vec(),
            })
        }
    }

    #[test]
    fn test_repeatmodeler_runs_with_mock() {
        let tmp = tempdir().unwrap();
        let fasta = tmp.path().join("genome.fa");
        std::fs::write(&fasta, ">x\nACGT").unwrap();

        let args = CliArgs {
            fasta_file: fasta,
            configure: Some(tmp.path().to_path_buf()),
            database: Some("mockdb".to_string()),
            rmo_threads: 1,
            rma_threads: 1,
            rma_only: false,
            verbose: false,
        };

        let runner = MockRunner;

        set_up_filesystem(args.clone()).unwrap();
        let out = run_repeatmodeler(args, &runner);
        assert!(out.is_ok());
    }
}
