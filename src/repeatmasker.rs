use std::path::{Path, PathBuf};
use std::process::Command;
use walkdir::WalkDir;

use crate::{
    command_runner::run_checked, configured_executable, CliArgs, CommandRunner, Error, ErrorKind,
    Result, DATA, REPEATMASKER_ENV,
};

pub fn run_repeatmasker(matches: CliArgs, runner: &dyn CommandRunner) -> Result<()> {
    // get the data path again
    let mut data_path = matches.configure.clone().unwrap();
    data_path.push(DATA);
    // and go into the repeatmasker dir
    data_path.push("repeatmasker");

    let mut repeatmodeler_path = matches.configure.clone().ok_or_else(|| {
        Error::new(ErrorKind::GenericCli(
            "A configuration directory is required".into(),
        ))
    })?;
    repeatmodeler_path.push(DATA);
    repeatmodeler_path.push("repeatmodeler");
    let full_consensi_path = find_consensus_library(&repeatmodeler_path)?;

    eprintln!("Data path: {:?}", data_path);

    let mut run_repeat_masker =
        Command::new(configured_executable(REPEATMASKER_ENV, "RepeatMasker"));
    run_repeat_masker
        // the number of threads
        .arg("-pa")
        .arg(matches.rma_threads.to_string())
        // the consensi.fa.classified file from the initial round of repeatmodeler
        .arg("-lib")
        .arg(full_consensi_path.to_string_lossy().to_string())
        // we want a gff
        .arg("-gff")
        // we want the alignment output
        .arg("-a")
        // and repeat densities
        .arg("-excln")
        .arg("-dir")
        .arg(&data_path)
        // and the genome file
        .arg(matches.fasta_file);

    run_checked(runner, &mut run_repeat_masker, "RepeatMasker")?;

    Ok(())
}

fn find_consensus_library(repeatmodeler_path: &Path) -> Result<PathBuf> {
    let mut matches = Vec::new();
    for entry in WalkDir::new(repeatmodeler_path).follow_links(false) {
        let entry = entry.map_err(|error| {
            Error::new(ErrorKind::GenericCli(format!(
                "Could not inspect RepeatModeler output: {error}"
            )))
        })?;
        if entry.file_type().is_file() && entry.file_name() == "consensi.fa.classified" {
            matches.push(entry.into_path());
        }
    }
    matches.sort();

    match matches.as_slice() {
        [path] => Ok(path.clone()),
        [] => Err(Error::new(ErrorKind::GenericCli(format!(
            "No consensi.fa.classified found under {}. Did RepeatModeler finish?",
            repeatmodeler_path.display()
        )))),
        paths => Err(Error::new(ErrorKind::GenericCli(format!(
            "Multiple consensi.fa.classified files found under {}: {}",
            repeatmodeler_path.display(),
            paths
                .iter()
                .map(|path| path.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        )))),
    }
}

#[cfg(test)]
mod tests {
    use crate::set_up_filesystem;

    use super::*;
    use std::fs::write;
    use std::os::unix::process::ExitStatusExt;
    use std::process::{Command, ExitStatus, Output};
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
    fn rejects_ambiguous_consensus_libraries() {
        let tmp = tempdir().unwrap();
        for directory in ["RM_1", "RM_2"] {
            let path = tmp.path().join(directory);
            std::fs::create_dir_all(&path).unwrap();
            write(path.join("consensi.fa.classified"), ">repeat\nACGT").unwrap();
        }

        let error = find_consensus_library(tmp.path()).unwrap_err();
        assert!(error
            .to_string()
            .contains("Multiple consensi.fa.classified files"));
    }

    #[test]
    fn test_repeatmasker_runs_with_mock() {
        let tmp = tempdir().unwrap();
        let fasta = tmp.path().join("genome.fa");
        write(&fasta, ">x\nACGT").unwrap();

        let consensi = tmp.path().join("consensi.fa.classified");
        write(&consensi, ">repeat\nACGT").unwrap();

        // simulate directory structure
        let data_dir = tmp.path().join("data").join("repeatmasker");
        std::fs::create_dir_all(&data_dir).unwrap();

        let args = CliArgs {
            fasta_file: fasta,
            configure: Some(tmp.path().to_path_buf()),
            database: Some("mockdb".to_string()),
            rmo_threads: 1,
            rma_threads: 1,
            rma_only: false,
            verbose: false,
        };

        set_up_filesystem(args.clone()).unwrap();

        // place the fake consensus in RepeatModeler output
        let consensi_path = tmp
            .path()
            .join("data")
            .join("repeatmodeler")
            .join("RM_test");
        std::fs::create_dir_all(&consensi_path).unwrap();
        let real_path = consensi_path.join("consensi.fa.classified");
        std::fs::copy(&consensi, &real_path).unwrap();

        // use the mock runner
        let runner = MockRunner;
        assert!(run_repeatmasker(args, &runner).is_ok());
    }
}
