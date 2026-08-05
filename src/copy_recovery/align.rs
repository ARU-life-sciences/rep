use std::path::Path;
use std::process::Command;
pub fn mafft_command(exe: impl AsRef<Path>, input: &Path, threads: usize) -> Command {
    let mut c = Command::new(exe.as_ref());
    c.args(["--auto", "--thread"])
        .arg(threads.to_string())
        .arg(input);
    c
}
