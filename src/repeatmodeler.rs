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
        .arg(&matches.fasta_file);

    run_checked(runner, &mut build_database, "BuildDatabase")?;

    let mut run_repeat_modeler =
        Command::new(configured_executable(REPEATMODELER_ENV, "RepeatModeler"));

    run_repeat_modeler
        .current_dir(data_path)
        .arg("-database")
        .arg(&database)
        .arg("-threads")
        .arg(matches.rmo_threads.to_string());
    if matches.ltr_struct {
        run_repeat_modeler.arg("-LTRStruct");
    }

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
            curation_target: None,
            configure: Some(tmp.path().to_path_buf()),
            database: Some("mockdb".to_string()),
            rmo_threads: 1,
            rma_threads: 1,
            rma_only: false,
            curation_only: false,
            verbose: false,
            ltr_struct: false,
            refdb: None,
        };

        let runner = MockRunner;

        set_up_filesystem(args.clone()).unwrap();
        let out = run_repeatmodeler(args, &runner);
        assert!(out.is_ok());
    }

    struct Recorder(std::cell::RefCell<Vec<Vec<String>>>);

    impl CommandRunner for Recorder {
        fn run(&self, cmd: &mut Command) -> Result<Output> {
            self.0.borrow_mut().push(
                std::iter::once(cmd.get_program())
                    .chain(cmd.get_args())
                    .map(|s| s.to_string_lossy().into_owned())
                    .collect(),
            );
            Ok(Output {
                status: ExitStatus::from_raw(0),
                stdout: Vec::new(),
                stderr: Vec::new(),
            })
        }
    }

    #[test]
    fn passes_ltr_struct_and_indexes_only_the_genome() {
        let tmp = tempdir().unwrap();
        let fasta = tmp.path().join("genome.fa");
        std::fs::write(&fasta, ">x\nACGT").unwrap();
        let args = CliArgs {
            fasta_file: fasta.clone(),
            curation_target: None,
            configure: Some(tmp.path().to_path_buf()),
            database: Some("db".to_string()),
            rmo_threads: 4,
            rma_threads: 1,
            rma_only: false,
            curation_only: false,
            verbose: false,
            ltr_struct: true,
            refdb: None,
        };
        set_up_filesystem(args.clone()).unwrap();
        let runner = Recorder(Default::default());
        run_repeatmodeler(args, &runner).unwrap();
        let calls = runner.0.into_inner();
        assert_eq!(&calls[0][1..], ["-name", "db", fasta.to_str().unwrap()]);
        assert_eq!(
            &calls[1][1..],
            ["-database", "db", "-threads", "4", "-LTRStruct"]
        );
    }
}
