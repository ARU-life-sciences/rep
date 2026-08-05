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
pub fn blastn_command(
    exe: impl AsRef<Path>,
    db: &Path,
    query: &Path,
    out: &Path,
    threads: usize,
) -> Command {
    let mut c = Command::new(exe.as_ref());
    c.args(["-db"]).arg(db).args(["-query"]).arg(query).args(["-outfmt", "6 qseqid sseqid qlen pident length mismatch gapopen qstart qend sstart send evalue bitscore", "-out"]).arg(out).args(["-num_threads"]).arg(threads.to_string());
    c
}
pub fn executable(name: &str) -> PathBuf {
    std::env::var_os(format!("REP_{}", name.to_uppercase()))
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(name))
}
