use std::process::Command;
use std::{fs, path::PathBuf};
use walkdir::WalkDir;

use crate::{CliArgs, Error, ErrorKind, Result, DATA};

pub fn run_repeatmasker(matches: CliArgs) -> Result<()> {
    // get the data path again
    let mut data_path = matches.configure.clone();
    data_path.push(DATA);
    // and go into the repeatmasker dir
    data_path.push("repeatmasker");

    // we need to find the consensi.fa.classified
    // inside the data directory
    let file_to_find = "consensi.fa.classified";
    let consensi_finder = || -> Option<PathBuf> {
        for entry in WalkDir::new(matches.configure.clone())
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

    let run_repeat_masker = Command::new("RepeatMasker")
        // run in the repeatmasker data subdirectory
        .current_dir(data_path)
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
        // and the genome file
        .arg(matches.fasta_file)
        .spawn()?;

    let _run_repeat_masker_out = run_repeat_masker.wait_with_output()?;

    Ok(())
}
