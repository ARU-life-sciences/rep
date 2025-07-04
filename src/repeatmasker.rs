use std::process::Command;
use std::{fs, path::PathBuf};
use walkdir::WalkDir;

use crate::{CliArgs, CommandRunner, Error, ErrorKind, Result, DATA};

pub fn run_repeatmasker(matches: CliArgs, runner: &dyn CommandRunner) -> Result<()> {
    // get the data path again
    let mut data_path = matches.configure.clone().unwrap();
    data_path.push(DATA);
    // and go into the repeatmasker dir
    data_path.push("repeatmasker");

    // we need to find the consensi.fa.classified
    // inside the data directory
    let file_to_find = "consensi.fa.classified";
    let consensi_finder = || -> Option<PathBuf> {
        for entry in WalkDir::new(matches.configure.clone().unwrap())
            .follow_links(true)
            .into_iter()
            .filter_map(|e| e.ok())
        {
            let f_name = entry.file_name().to_string_lossy();

            if f_name == file_to_find {
                return Some(fs::canonicalize(entry.path()).unwrap());
            }
        }
        None
    };

    let full_consensi_path =
        match consensi_finder() {
            Some(p) => p,
            None => return Err(Error::new(ErrorKind::GenericCli(
                "No consensi.fa.classified found in the data directory. Did you run RepeatModeler?"
                    .to_string(),
            ))),
        };

    eprintln!("Data path: {:?}", data_path);

    let mut run_repeat_masker =
        Command::new("/software/team301/repeat-annotation/RepeatMasker/RepeatMasker");
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

    let output = runner.run(&mut run_repeat_masker)?;
    if !output.status.success() {
        return Err(Error::new(ErrorKind::GenericCli(format!(
            "RepeatMasker failed: {}",
            String::from_utf8_lossy(&output.stderr)
        ))));
    }

    Ok(())
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
            curation_only: None,
            curation_rmdl_library: None,
            verbose: false,
        };

        set_up_filesystem(args.clone()).unwrap();

        // place the fake consensi file somewhere discoverable
        let consensi_path = tmp.path().join("dummy_consensi_dir");
        std::fs::create_dir_all(&consensi_path).unwrap();
        let real_path = consensi_path.join("consensi.fa.classified");
        std::fs::copy(&consensi, &real_path).unwrap();

        // use the mock runner
        let runner = MockRunner;
        assert!(run_repeatmasker(args, &runner).is_ok());
    }
}
