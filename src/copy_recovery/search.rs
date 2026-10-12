use std::path::{Path, PathBuf};
use std::process::Command;

pub fn makeblastdb_command(exe: impl AsRef<Path>, genome: &Path, out: &Path) -> Command {
    let mut c = Command::new(exe.as_ref());
    c.args(["-in"])
        .arg(genome)
        .args(["-dbtype", "nucl", "-out"])
        .arg(out)
        .arg("-parse_seqids");
    c
}
/// Build the BLASTN genome search. `task` should be a sensitive task
/// (`dc-megablast` or `blastn`): the megablast default misses diverged copies.
pub fn blastn_command(
    exe: impl AsRef<Path>,
    db: &Path,
    query: &Path,
    out: &Path,
    threads: usize,
    task: &str,
    evalue: f64,
) -> Command {
    let mut c = Command::new(exe.as_ref());
    c.args(["-task", task, "-db"])
        .arg(db)
        .args(["-query"])
        .arg(query)
        .args([
            "-outfmt",
            "6 qseqid sseqid qlen pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-evalue",
        ])
        .arg(evalue.to_string())
        .args(["-out"])
        .arg(out)
        .args(["-num_threads"])
        .arg(threads.to_string());
    c
}
pub fn executable(name: &str) -> PathBuf {
    std::env::var_os(format!("REP_{}", name.to_uppercase()))
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(name))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blastn_uses_sensitive_task_and_evalue() {
        let c = blastn_command(
            "blastn",
            Path::new("db"),
            Path::new("q.fa"),
            Path::new("o.tsv"),
            4,
            "dc-megablast",
            1e-10,
        );
        let args = c
            .get_args()
            .map(|a| a.to_string_lossy().into_owned())
            .collect::<Vec<_>>();
        assert_eq!(&args[..2], ["-task", "dc-megablast"]);
        let e = args.iter().position(|a| a == "-evalue").unwrap();
        assert_eq!(args[e + 1].parse::<f64>().unwrap(), 1e-10);
    }
}
