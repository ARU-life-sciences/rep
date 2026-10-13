//! Superfamily profile HMMs built from the labelled reference proteins.
//!
//! Pfam detects some TE superfamilies poorly (Crypton, Penelope, DIRS,
//! Maverick, P). Every labelled protein in the build is clustered with
//! MMseqs2; each cluster (split by label) with enough members is aligned
//! with MAFFT and turned into a profile with `hmmbuild`. Model names are
//! `hmm_<cluster>#<Class/Superfamily>`, so a hit carries its class.

use std::{
    collections::BTreeMap,
    fs,
    io::Write,
    path::{Path, PathBuf},
    process::Command,
    sync::Mutex,
};

use bio::io::fasta;
use serde::{Deserialize, Serialize};

use super::{cli_error, BuildInfo, RefDb};
use crate::{
    classify::{is_non_te, split_label},
    command_runner::run_checked,
    copy_recovery::executable,
    CommandRunner, Result,
};

#[derive(Clone, Debug, Serialize)]
pub struct HmmBuildOptions {
    pub identity: f64,
    pub coverage: f64,
    pub min_members: usize,
    /// Members aligned per model (longest first).
    pub max_members: usize,
    pub threads: usize,
}

impl Default for HmmBuildOptions {
    fn default() -> Self {
        Self {
            identity: 0.3,
            coverage: 0.5,
            min_members: 5,
            max_members: 60,
            threads: 1,
        }
    }
}

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct HmmBuildInfo {
    pub proteins: usize,
    pub clusters: usize,
    pub models: usize,
    pub labels: usize,
    pub failed: usize,
    /// `build.json` fingerprint of the protein library the models came from.
    pub source_fingerprint: String,
}

impl RefDb {
    pub fn hmm_dir(&self) -> PathBuf {
        self.root.join("hmms")
    }
    pub fn te_hmms(&self) -> PathBuf {
        self.hmm_dir().join("te_hmms.hmm")
    }
}

/// Build `<db>/hmms/te_hmms.hmm` (pressed) from `build/te_prot.fa`.
pub fn build_te_hmms(
    db: &RefDb,
    options: &HmmBuildOptions,
    runner: &(dyn CommandRunner + Sync),
) -> Result<HmmBuildInfo> {
    build_te_hmms_from(&db.protein_library(), &db.hmm_dir(), options, runner)
}

/// As [`build_te_hmms`] for any labelled protein FASTA (`id#label` IDs);
/// used for held-out benchmarks.
pub fn build_te_hmms_from(
    proteins: &Path,
    out_dir: &Path,
    options: &HmmBuildOptions,
    runner: &(dyn CommandRunner + Sync),
) -> Result<HmmBuildInfo> {
    let work = out_dir.join("work");
    if work.exists() {
        fs::remove_dir_all(&work)?;
    }
    fs::create_dir_all(work.join("models"))?;

    // usable labelled proteins, renamed p<index>
    let mut labels = Vec::new();
    let mut seqs = Vec::new();
    let mut fa = String::new();
    for record in fasta::Reader::from_file(proteins)?.records() {
        let record = record?;
        let (_, label) = split_label(record.id());
        let Some(label) =
            label.filter(|l| l.contains('/') && !l.contains("Unknown") && !is_non_te(l))
        else {
            continue;
        };
        if record.seq().len() < 50 {
            continue;
        }
        fa.push_str(&format!(
            ">p{}\n{}\n",
            labels.len(),
            String::from_utf8_lossy(record.seq())
        ));
        labels.push(label);
        seqs.push(record.seq().to_vec());
    }
    let input = work.join("proteins.fa");
    fs::write(&input, fa)?;
    // a refdb protein library has its build record beside it
    let source_fingerprint = fs::read(proteins.with_file_name("build.json"))
        .ok()
        .and_then(|b| serde_json::from_slice::<BuildInfo>(&b).ok())
        .map(|b| b.fingerprint)
        .unwrap_or_default();
    let mut info = HmmBuildInfo {
        proteins: labels.len(),
        source_fingerprint,
        ..Default::default()
    };

    let mut cluster = Command::new(executable("mmseqs"));
    cluster
        .arg("easy-cluster")
        .arg(&input)
        .arg(work.join("clu"))
        .arg(work.join("tmp"))
        .args([
            "--min-seq-id",
            &options.identity.to_string(),
            "-c",
            &options.coverage.to_string(),
            "--cov-mode",
            "0",
            "-v",
            "1",
            "--threads",
            &options.threads.to_string(),
        ]);
    run_checked(runner, &mut cluster, "mmseqs easy-cluster")?;
    let index = |id: &str| id.trim_start_matches('p').parse::<usize>().ok();
    let mut groups: BTreeMap<(usize, String), Vec<usize>> = BTreeMap::new();
    for line in fs::read_to_string(work.join("clu_cluster.tsv"))?.lines() {
        if let Some((rep, member)) = line.split_once('\t') {
            if let (Some(r), Some(m)) = (index(rep), index(member)) {
                groups.entry((r, labels[m].clone())).or_default().push(m);
            }
        }
    }
    info.clusters = groups.len();
    let mut jobs = groups
        .into_iter()
        .filter(|(_, members)| members.len() >= options.min_members)
        .map(|((rep, label), mut members)| {
            members.sort_by_key(|&m| (std::cmp::Reverse(seqs[m].len()), m));
            members.truncate(options.max_members);
            (rep, label, members)
        })
        .collect::<Vec<_>>();
    jobs.sort_by_key(|(rep, label, _)| (label.clone(), *rep));

    // align and build models in parallel; each worker uses one CPU
    let results = Mutex::new(Vec::new());
    let next = Mutex::new(0usize);
    std::thread::scope(|scope| {
        for _ in 0..options.threads.max(1) {
            scope.spawn(|| loop {
                let i = {
                    let mut n = next.lock().expect("lock");
                    let i = *n;
                    *n += 1;
                    i
                };
                let Some((rep, label, members)) = jobs.get(i) else {
                    break;
                };
                let name = format!("hmm_{rep}#{label}");
                let built = build_one(&work, &name, members, &seqs, runner);
                results
                    .lock()
                    .expect("lock")
                    .push((i, label.clone(), built));
            });
        }
    });
    let mut results = results.into_inner().expect("lock");
    results.sort_by_key(|(i, _, _)| *i);

    let tmp = out_dir.join(format!(".te_hmms.hmm.tmp-{}", std::process::id()));
    let mut out = fs::File::create(&tmp)?;
    let mut table = String::from("model\tlabel\tmembers\n");
    let mut built_labels = std::collections::BTreeSet::new();
    for (i, label, built) in results {
        match built {
            Ok(hmm) => {
                out.write_all(hmm.as_bytes())?;
                let (rep, _, members) = &jobs[i];
                table.push_str(&format!("hmm_{rep}#{label}\t{label}\t{}\n", members.len()));
                built_labels.insert(label);
                info.models += 1;
            }
            Err(e) => {
                eprintln!("model {i} ({label}) failed: {e}");
                info.failed += 1;
            }
        }
    }
    out.flush()?;
    drop(out);
    if info.models == 0 {
        return Err(cli_error("no TE profile HMMs could be built"));
    }
    info.labels = built_labels.len();
    let target = out_dir.join("te_hmms.hmm");
    for suffix in ["", ".h3f", ".h3i", ".h3m", ".h3p"] {
        fs::remove_file(format!("{}{suffix}", target.display())).ok();
    }
    fs::rename(&tmp, &target)?;
    let mut press = Command::new(executable("hmmpress"));
    press.arg("-f").arg(&target);
    run_checked(runner, &mut press, "hmmpress")?;
    fs::write(out_dir.join("te_hmms.tsv"), table)?;
    fs::write(
        out_dir.join("te_hmms.json"),
        serde_json::to_vec_pretty(&info).map_err(|e| cli_error(e.to_string()))?,
    )?;
    fs::remove_dir_all(&work).ok();
    Ok(info)
}

fn build_one(
    work: &Path,
    name: &str,
    members: &[usize],
    seqs: &[Vec<u8>],
    runner: &dyn CommandRunner,
) -> Result<String> {
    let tag = name.replace(['#', '/'], "_");
    let input = work.join("models").join(format!("{tag}.fa"));
    let mut fa = String::new();
    for &m in members {
        fa.push_str(&format!(">p{m}\n{}\n", String::from_utf8_lossy(&seqs[m])));
    }
    fs::write(&input, fa)?;
    let mut mafft = Command::new(executable("mafft"));
    mafft
        .args(["--auto", "--anysymbol", "--quiet", "--thread", "1"])
        .arg(&input);
    let aligned = run_checked(runner, &mut mafft, "mafft")?;
    let aln = work.join("models").join(format!("{tag}.aln"));
    fs::write(&aln, aligned.stdout)?;
    let hmm = work.join("models").join(format!("{tag}.hmm"));
    let mut build = Command::new(executable("hmmbuild"));
    build
        .args(["--amino", "--cpu", "1", "-o", "/dev/null", "-n", name])
        .arg(&hmm)
        .arg(&aln);
    run_checked(runner, &mut build, "hmmbuild")?;
    let text = fs::read_to_string(&hmm)?;
    fs::remove_file(&input).ok();
    fs::remove_file(&aln).ok();
    fs::remove_file(&hmm).ok();
    Ok(text)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::process::ExitStatusExt;
    use std::process::{ExitStatus, Output};

    /// Fakes mmseqs (two clusters), mafft (echo input) and hmmbuild (tiny model).
    struct Fake;

    impl CommandRunner for Fake {
        fn run(&self, cmd: &mut Command) -> Result<Output> {
            let prog = cmd.get_program().to_string_lossy().into_owned();
            let args = cmd
                .get_args()
                .map(|a| a.to_string_lossy().into_owned())
                .collect::<Vec<_>>();
            let mut stdout = Vec::new();
            if prog.ends_with("mmseqs") {
                // p0..p5 Gypsy in one cluster with a mislabelled p6; p7..p8 too small
                let rows =
                    (0..7).map(|m| format!("p0\tp{m}\n")).collect::<String>() + "p7\tp7\np7\tp8\n";
                fs::write(format!("{}_cluster.tsv", args[2]), rows)?;
            } else if prog.ends_with("mafft") {
                stdout = fs::read(args.last().unwrap())?;
            } else if prog.ends_with("hmmbuild") {
                let name = &args[args.iter().position(|a| a == "-n").unwrap() + 1];
                fs::write(
                    &args[args.len() - 2],
                    format!("HMMER3/f\nNAME  {name}\n//\n"),
                )?;
            }
            Ok(Output {
                status: ExitStatus::from_raw(0),
                stdout,
                stderr: vec![],
            })
        }
    }

    #[test]
    fn builds_labelled_models_from_clusters() {
        let dir = tempfile::tempdir().unwrap();
        let prot = dir.path().join("p.fa");
        let seq = "M".repeat(80);
        let mut fa = String::new();
        for i in 0..6 {
            fa.push_str(&format!(">g{i}#LTR/Gypsy\n{seq}\n"));
        }
        fa.push_str(&format!(">x#DNA/hAT\n{seq}\n>y#DNA/hAT\n{seq}\n>z#DNA/hAT\n{seq}\n>u#Unknown\n{seq}\n>short#LTR/Gypsy\nMK\n"));
        fs::write(&prot, fa).unwrap();
        let out = dir.path().join("hmms");
        let info = build_te_hmms_from(&prot, &out, &HmmBuildOptions::default(), &Fake).unwrap();
        // Unknown and too-short proteins are excluded
        assert_eq!(info.proteins, 9);
        // the Gypsy cluster gives one model; hAT members are split off and too few
        assert_eq!((info.models, info.labels, info.failed), (1, 1, 0));
        let hmm = fs::read_to_string(out.join("te_hmms.hmm")).unwrap();
        assert_eq!(hmm, "HMMER3/f\nNAME  hmm_0#LTR/Gypsy\n//\n");
        assert!(fs::read_to_string(out.join("te_hmms.tsv"))
            .unwrap()
            .contains("hmm_0#LTR/Gypsy\tLTR/Gypsy\t6"));
        assert!(!out.join("work").exists());
        assert_eq!(info.source_fingerprint, "");
        fs::write(dir.path().join("build.json"), r#"{"fingerprint":"abc"}"#).unwrap();
        let info = build_te_hmms_from(&prot, &out, &HmmBuildOptions::default(), &Fake).unwrap();
        assert_eq!(info.source_fingerprint, "abc");
    }
}
