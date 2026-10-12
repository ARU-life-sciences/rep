//! Download and register public source releases.
//!
//! Downloads go through the [`CommandRunner`] (`curl`, `tar`, `gzip`,
//! `md5sum`) so the logic is testable without network access. Releases that
//! are already registered with the same checksum are skipped, which makes
//! `fetch` an idempotent update.

use std::{fs, path::Path, process::Command};

use serde::Deserialize;

use super::{
    cli_error, now, sha256_file, RefDb, SourceRelease, DFAM_CURATED, DFAM_PEPTIDES, REPBASE,
};
use crate::{command_runner::run_checked, CommandRunner, Result};

pub const REPBASE_MANIFEST_URL: &str = "https://www.dfam.org/repbase_archive/releases.json";
pub const REPBASE_ARCHIVE_URL: &str = "https://www.dfam.org/repbase_archive";
pub const DFAM_CURRENT_URL: &str = "https://www.dfam.org/releases/current";
pub const CC0: &str = "CC0-1.0";

#[derive(Clone, Debug)]
pub struct FetchOptions {
    pub repbase: bool,
    /// A specific Repbase release such as `RepBase31.09`; latest when `None`.
    pub repbase_release: Option<String>,
    pub dfam: bool,
}

impl Default for FetchOptions {
    fn default() -> Self {
        Self {
            repbase: true,
            repbase_release: None,
            dfam: true,
        }
    }
}

#[derive(Debug, Deserialize)]
struct RepbaseManifest {
    releases: Vec<RepbaseRelease>,
}

#[derive(Debug, Deserialize)]
struct RepbaseRelease {
    release: String,
    files: RepbaseFiles,
}

#[derive(Debug, Deserialize)]
struct RepbaseFiles {
    #[serde(default)]
    embl: Vec<RepbaseFile>,
}

#[derive(Debug, Deserialize)]
struct RepbaseFile {
    name: String,
    path: String,
    sha256: String,
}

/// Fetch the selected sources; returns the releases newly registered.
pub fn fetch(
    db: &RefDb,
    options: &FetchOptions,
    runner: &dyn CommandRunner,
) -> Result<Vec<SourceRelease>> {
    fs::create_dir_all(db.root.join("cache"))?;
    let mut added = Vec::new();
    if options.repbase {
        added.extend(fetch_repbase(
            db,
            options.repbase_release.as_deref(),
            runner,
        )?);
    }
    if options.dfam {
        added.extend(fetch_dfam(db, runner)?);
    }
    Ok(added)
}

fn download(runner: &dyn CommandRunner, url: &str, out: &Path) -> Result<()> {
    let mut cmd = Command::new("curl");
    cmd.args(["-fsSL", "--retry", "3", "-o"]).arg(out).arg(url);
    run_checked(runner, &mut cmd, &format!("download of {url}"))?;
    Ok(())
}

fn fetch_repbase(
    db: &RefDb,
    wanted: Option<&str>,
    runner: &dyn CommandRunner,
) -> Result<Option<SourceRelease>> {
    let manifest_path = db.root.join("cache").join("repbase_releases.json");
    download(runner, REPBASE_MANIFEST_URL, &manifest_path)?;
    let releases: RepbaseManifest = serde_json::from_slice(&fs::read(&manifest_path)?)
        .map_err(|e| cli_error(format!("invalid Repbase release manifest: {e}")))?;
    // the manifest lists releases newest first
    let release = match wanted {
        Some(name) => releases
            .releases
            .iter()
            .find(|r| r.release == name)
            .ok_or_else(|| cli_error(format!("Repbase release {name} not found")))?,
        None => releases
            .releases
            .first()
            .ok_or_else(|| cli_error("Repbase release manifest is empty"))?,
    };
    let file = release
        .files
        .embl
        .first()
        .ok_or_else(|| cli_error(format!("{} has no EMBL file", release.release)))?;
    let version = release.release.trim_start_matches("RepBase").to_owned();
    if already_registered(db, REPBASE, &version, &file.sha256)? {
        eprintln!("Repbase {version} is already registered");
        return Ok(None);
    }
    let dir = db.root.join("sources").join(REPBASE).join(&version);
    fs::create_dir_all(&dir)?;
    let archive = dir.join(&file.name);
    let url = format!("{REPBASE_ARCHIVE_URL}/{}", file.path);
    eprintln!("Downloading Repbase {version} from {url}");
    download(runner, &url, &archive)?;
    let actual = sha256_file(&archive)?;
    if actual != file.sha256 {
        return Err(cli_error(format!(
            "SHA-256 mismatch for {}: expected {}, got {actual}",
            archive.display(),
            file.sha256
        )));
    }
    register_source(db, REPBASE, &version, &archive, Some(url), runner).map(Some)
}

fn fetch_dfam(db: &RefDb, runner: &dyn CommandRunner) -> Result<Vec<SourceRelease>> {
    let notes = db.root.join("cache").join("dfam_relnotes.txt");
    download(runner, &format!("{DFAM_CURRENT_URL}/relnotes.txt"), &notes)?;
    let version = dfam_version(&fs::read_to_string(&notes)?)
        .ok_or_else(|| cli_error("could not find the Dfam release number in relnotes.txt"))?;
    let mut added = Vec::new();
    for (name, file) in [
        (DFAM_CURATED, "Dfam-RepeatMasker.lib.gz"),
        (DFAM_PEPTIDES, "RepeatPeps.lib.gz"),
    ] {
        let dir = db.root.join("sources").join(name).join(&version);
        fs::create_dir_all(&dir)?;
        let url = format!("{DFAM_CURRENT_URL}/families/{file}");
        let md5_path = dir.join(format!("{file}.md5"));
        download(runner, &format!("{url}.md5"), &md5_path)?;
        let expected = fs::read_to_string(&md5_path)?
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_owned();
        if db
            .load_manifest()?
            .sources
            .iter()
            .any(|s| s.name == name && s.version == version)
        {
            eprintln!("{name} {version} is already registered");
            continue;
        }
        let archive = dir.join(file);
        eprintln!("Downloading {name} {version} from {url}");
        download(runner, &url, &archive)?;
        let mut md5 = Command::new("md5sum");
        md5.arg(&archive);
        let out = run_checked(runner, &mut md5, "md5sum")?;
        let actual = String::from_utf8_lossy(&out.stdout)
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_owned();
        if actual != expected {
            return Err(cli_error(format!(
                "MD5 mismatch for {}: expected {expected}, got {actual}",
                archive.display()
            )));
        }
        added.push(register_source(
            db,
            name,
            &version,
            &archive,
            Some(url),
            runner,
        )?);
    }
    Ok(added)
}

/// `RELEASE 4.0` in the Dfam release notes.
fn dfam_version(notes: &str) -> Option<String> {
    notes.lines().find_map(|line| {
        let mut words = line.split_whitespace();
        (words.next() == Some("RELEASE"))
            .then(|| words.next().map(str::to_owned))
            .flatten()
    })
}

fn already_registered(db: &RefDb, name: &str, version: &str, sha256: &str) -> Result<bool> {
    Ok(db
        .load_manifest()?
        .sources
        .iter()
        .any(|s| s.name == name && s.version == version && s.sha256.iter().any(|h| h == sha256)))
}

/// Register a local file as a source release: `.tar.gz` archives are unpacked
/// and `.gz` files decompressed inside the database. Repbase expects an EMBL
/// archive, file or directory; Dfam sources expect labelled FASTA.
pub fn register_source(
    db: &RefDb,
    name: &str,
    version: &str,
    file: &Path,
    url: Option<String>,
    runner: &dyn CommandRunner,
) -> Result<SourceRelease> {
    if ![REPBASE, DFAM_CURATED, DFAM_PEPTIDES].contains(&name) {
        return Err(cli_error(format!(
            "unknown source {name}; expected {REPBASE}, {DFAM_CURATED} or {DFAM_PEPTIDES}"
        )));
    }
    let dir = db.root.join("sources").join(name).join(version);
    fs::create_dir_all(&dir)?;
    let file_name = file
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let sha256 = if file.is_file() {
        vec![sha256_file(file)?]
    } else {
        Vec::new()
    };
    let unpacked = if file_name.ends_with(".tar.gz") || file_name.ends_with(".tgz") {
        let target = dir.join("unpacked");
        fs::create_dir_all(&target)?;
        let mut tar = Command::new("tar");
        tar.arg("-xzf").arg(file).arg("-C").arg(&target);
        run_checked(runner, &mut tar, "tar")?;
        target
    } else if let Some(stem) = file_name.strip_suffix(".gz") {
        let mut gzip = Command::new("gzip");
        gzip.arg("-dc").arg(file);
        let out = run_checked(runner, &mut gzip, "gzip")?;
        let target = dir.join(stem);
        fs::write(&target, out.stdout)?;
        target
    } else if file.starts_with(&dir) {
        file.to_path_buf()
    } else if file.is_dir() {
        return Err(cli_error(
            "register a directory by archiving it, or place it under the database sources/ directory",
        ));
    } else {
        let target = dir.join(&file_name);
        fs::copy(file, &target)?;
        target
    };
    let relative = unpacked
        .strip_prefix(&db.root)
        .map(|p| p.to_string_lossy().into_owned())
        .map_err(|_| cli_error("source files must live inside the database"))?;
    let release = SourceRelease {
        name: name.into(),
        version: version.into(),
        url,
        files: vec![relative],
        sha256,
        licence: CC0.into(),
        registered_at: now(),
    };
    let mut manifest = db.load_manifest()?;
    manifest
        .sources
        .retain(|s| !(s.name == name && s.version == version));
    manifest.sources.push(release.clone());
    db.save_manifest(&manifest)?;
    Ok(release)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::process::ExitStatusExt;
    use std::process::{ExitStatus, Output};

    /// Serves canned files for `curl -o <out> <url>` and runs other tools for real.
    struct FakeNet {
        files: Vec<(&'static str, Vec<u8>)>,
    }

    impl CommandRunner for FakeNet {
        fn run(&self, cmd: &mut Command) -> Result<Output> {
            if cmd.get_program() != "curl" {
                return crate::RealCommandRunner.run(cmd);
            }
            let args = cmd
                .get_args()
                .map(|a| a.to_string_lossy().into_owned())
                .collect::<Vec<_>>();
            let out = &args[args.iter().position(|a| a == "-o").unwrap() + 1];
            let url = args.last().unwrap();
            let body = self
                .files
                .iter()
                .find(|(u, _)| url.ends_with(u))
                .unwrap_or_else(|| panic!("unexpected download {url}"));
            fs::write(out, &body.1)?;
            Ok(Output {
                status: ExitStatus::from_raw(0),
                stdout: Vec::new(),
                stderr: Vec::new(),
            })
        }
    }

    fn gzip(data: &[u8]) -> Vec<u8> {
        let mut child = Command::new("gzip")
            .arg("-c")
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .spawn()
            .unwrap();
        use std::io::Write;
        child.stdin.take().unwrap().write_all(data).unwrap();
        child.wait_with_output().unwrap().stdout
    }

    #[test]
    fn parses_dfam_release_number() {
        assert_eq!(
            dfam_version("  header\n   RELEASE 4.0\n").as_deref(),
            Some("4.0")
        );
        assert_eq!(dfam_version("nothing"), None);
    }

    #[test]
    fn fetches_verifies_registers_and_skips_existing_releases() {
        let dir = tempfile::tempdir().unwrap();
        let db = RefDb::new(dir.path().join("db"));
        // a Repbase EMBL archive
        let src = dir.path().join("src");
        fs::create_dir_all(src.join("RB")).unwrap();
        fs::write(
            src.join("RB/x.ref"),
            "ID   A_X   repbase;\nKW   hAT; DNA transposon;\nSQ   Sequence\n     acgt 4\n//\n",
        )
        .unwrap();
        let archive = dir.path().join("RB.embl.tar.gz");
        assert!(Command::new("tar")
            .arg("-czf")
            .arg(&archive)
            .arg("-C")
            .arg(&src)
            .arg("RB")
            .status()
            .unwrap()
            .success());
        let archive_bytes = fs::read(&archive).unwrap();
        let sha = super::super::sha256_hex(&archive_bytes);
        let manifest = format!(
            r#"{{"releases":[{{"release":"RepBase31.09","date":"x","files":{{"embl":[{{"name":"RB.embl.tar.gz","path":"RepBase31.09/RB.embl.tar.gz","sha256":"{sha}"}}]}}}}]}}"#
        );
        let lib = gzip(b">F1#DNA/hAT @Clade\nAAAACCCCGGGG\n");
        let peps = gzip(b">P1#LTR/Gypsy\nMKV\n");
        let md5 = |d: &[u8]| {
            let p = dir.path().join("tmp.bin");
            fs::write(&p, d).unwrap();
            let o = Command::new("md5sum").arg(&p).output().unwrap();
            String::from_utf8(o.stdout).unwrap()
        };
        let net = FakeNet {
            files: vec![
                ("releases.json", manifest.into_bytes()),
                ("RB.embl.tar.gz", archive_bytes),
                ("relnotes.txt", b"   RELEASE 4.0\n".to_vec()),
                ("Dfam-RepeatMasker.lib.gz.md5", md5(&lib).into_bytes()),
                ("Dfam-RepeatMasker.lib.gz", lib),
                ("RepeatPeps.lib.gz.md5", md5(&peps).into_bytes()),
                ("RepeatPeps.lib.gz", peps),
            ],
        };
        let added = fetch(&db, &FetchOptions::default(), &net).unwrap();
        assert_eq!(
            added
                .iter()
                .map(|r| (r.name.as_str(), r.version.as_str()))
                .collect::<Vec<_>>(),
            vec![
                (REPBASE, "31.09"),
                (DFAM_CURATED, "4.0"),
                (DFAM_PEPTIDES, "4.0")
            ]
        );
        assert!(fetch(&db, &FetchOptions::default(), &net)
            .unwrap()
            .is_empty());

        let info = super::super::build(&db).unwrap();
        assert_eq!((info.nucleotide_records, info.protein_records), (2, 1));
        let nucl = fs::read_to_string(db.nucleotide_library()).unwrap();
        assert!(nucl.contains(">A_X#DNA/hAT source=repbase@31.09\n"));
        assert!(nucl.contains(">F1#DNA/hAT source=dfam-curated@4.0 [Clade]\n"));
    }

    #[test]
    fn rejects_checksum_mismatches() {
        let dir = tempfile::tempdir().unwrap();
        let db = RefDb::new(dir.path());
        let manifest = r#"{"releases":[{"release":"RepBase1.0","files":{"embl":[{"name":"a.tar.gz","path":"a.tar.gz","sha256":"00"}]}}]}"#;
        let net = FakeNet {
            files: vec![
                ("releases.json", manifest.as_bytes().to_vec()),
                ("a.tar.gz", b"not the archive".to_vec()),
            ],
        };
        let options = FetchOptions {
            dfam: false,
            ..Default::default()
        };
        let err = fetch(&db, &options, &net).unwrap_err();
        assert!(err.to_string().contains("SHA-256 mismatch"));
    }
}
