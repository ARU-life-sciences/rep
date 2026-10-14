use clap::{arg, command, value_parser, ArgAction};

use crate::{error::Result, Error, ErrorKind};
use std::path::{Path, PathBuf};

// a struct to contain all the CliArgs
// at the moment we only want the path
// to the fasta file
#[derive(Debug, Clone)]
pub struct CliArgs {
    // path to the primary FASTA file
    pub fasta_file: PathBuf,
    // optional separate genome searched during consensus curation
    pub curation_target: Option<PathBuf>,
    // whether to configure file system
    pub configure: Option<PathBuf>,
    // the name of the database for BuildDatabase
    pub database: Option<String>,
    // repeat modeler threads
    pub rmo_threads: u8,
    // repeat masker threads
    pub rma_threads: u8,
    // run repeat masker only
    pub rma_only: bool,
    // run consensus curation from existing RepeatModeler output
    pub curation_only: bool,
    // verbose flag
    pub verbose: bool,
    // run RepeatModeler's LTR structural discovery (-LTRStruct)
    pub ltr_struct: bool,
    // reference database used to classify families before masking
    pub refdb: Option<PathBuf>,
    // Pfam-A.hmm for protein-domain architecture during classification
    pub pfam: Option<PathBuf>,
}

impl CliArgs {
    pub fn configure_path(&self) -> Result<&Path> {
        self.configure.as_deref().ok_or_else(|| {
            Error::new(ErrorKind::GenericCli(
                "A configuration directory is required".into(),
            ))
        })
    }

    pub fn database_name(&self) -> Result<&str> {
        self.database.as_deref().ok_or_else(|| {
            Error::new(ErrorKind::GenericCli(
                "A database name is required to run RepeatModeler".into(),
            ))
        })
    }
}

pub fn parse_args() -> Result<CliArgs> {
    let matches = command!()
        .next_line_help(true)
        // not optional
        .arg(
            arg!(<FASTA> "Input file in FASTA format.")
                .value_parser(value_parser!(PathBuf)
                )
        )
        .arg(
            arg!(--curation_target <CURATION_FASTA> "Optional separate FASTA searched for consensus copies.")
                .value_parser(value_parser!(PathBuf)),
        )
        // not optional
        .arg(
            arg!(-c --configure <CONFIG_PATH> "Configure the file system - and create the required directories.")
                .required(true)
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            arg!(-d --database <DATABASE_NAME> "Name of the database, when building using `BuildDatabase`.")
                // not required unless you only want to run RepeatMasker
                // or I think, when you want run the curation pipeline
                .required_unless_present_any(["rma_only", "curation_only"])
                .value_parser(value_parser!(String)),
        )
        .arg(
            arg!(--rmo_threads <RMO_THREADS> "Number of threads to use for RepeatModeler.")
                .default_value("8")
                .value_parser(value_parser!(u8)),

        )
        .arg(
            arg!(--rma_threads <RMA_THREADS> "Number of threads to use for RepeatMasker.")
                .default_value("8")
                .value_parser(value_parser!(u8)),

        )
        .arg(
            arg!(--rma_only "Run RepeatMasker only. Skip RepeatModeler and curation.")
                .conflicts_with("curation_only")
                .action(ArgAction::SetTrue)
        )
        .arg(
            arg!(--curation_only "Run curation from existing RepeatModeler output.")
                .action(ArgAction::SetTrue)
        )
        .arg(
            arg!(--ltr_struct "Run RepeatModeler with -LTRStruct (LTR structural discovery).")
                .action(ArgAction::SetTrue),
        )
        .arg(
            arg!(--refdb <REFDB> "Built `rep refdb` database: classify families and mask with the reclassified library.")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            arg!(--pfam <PFAM> "Pfam-A.hmm used for protein-domain architecture during classification (with --refdb).")
                .requires("refdb")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(arg!(--verbose "Print extra debug information").action(ArgAction::SetTrue))
        .get_matches();

    // parse the arguments out
    let fasta = matches
        .get_one::<PathBuf>("FASTA")
        .cloned()
        .expect("errored by clap");

    if !fasta.is_file() {
        return Err(Error::new(ErrorKind::GenericCli(format!(
            "FASTA path is not a regular file: {}",
            fasta.display()
        ))));
    }

    let curation_target = matches.get_one::<PathBuf>("curation_target").cloned();
    if let Some(path) = &curation_target {
        if !path.is_file() {
            return Err(Error::new(ErrorKind::GenericCli(format!(
                "Curation target is not a regular file: {}",
                path.display()
            ))));
        }
    }

    let configure = matches.get_one::<PathBuf>("configure").cloned();

    let database = matches.get_one::<String>("database").cloned();

    let rmo_threads = matches
        .get_one::<u8>("rmo_threads")
        .cloned()
        .expect("errored by clap");

    let rma_threads = matches
        .get_one::<u8>("rma_threads")
        .cloned()
        .expect("errored by clap");

    let rma_only = matches.get_flag("rma_only");

    let curation_only = matches.get_flag("curation_only");

    let verbose = matches.get_flag("verbose");

    let ltr_struct = matches.get_flag("ltr_struct");

    let refdb = matches.get_one::<PathBuf>("refdb").cloned();
    if let Some(db) = &refdb {
        if !db.join("build").join("te_nucl.fa").is_file() {
            return Err(Error::new(ErrorKind::GenericCli(format!(
                "{} is not a built reference database (run `rep refdb build`)",
                db.display()
            ))));
        }
    }

    // collect the arguments
    Ok(CliArgs {
        fasta_file: fasta,
        curation_target,
        configure,
        database,
        rmo_threads,
        rma_threads,
        rma_only,
        curation_only,
        verbose,
        ltr_struct,
        refdb,
        pfam: matches.get_one::<PathBuf>("pfam").cloned(),
    })
}

fn recover_copies_command() -> clap::Command {
    use clap::{Arg, Command};
    Command::new("recover-copies")
        .about("Recover, extract and align genomic copies of candidate repeat families.")
        .arg(
            Arg::new("queries")
                .long("queries")
                .required(true)
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("genome")
                .long("genome")
                .required(true)
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("out")
                .long("out")
                .required(true)
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("flank-length")
                .long("flank-length")
                .default_value("2000")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("hits-per-family")
                .long("hits-per-family")
                .default_value("20")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("max-fragment-gap")
                .long("max-fragment-gap")
                .default_value("10000")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("orientation-threshold")
                .long("orientation-threshold")
                .default_value("0.8")
                .value_parser(value_parser!(f64)),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .default_value("1")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("search-backend")
                .long("search-backend")
                .default_value("blastn"),
        )
        .arg(
            Arg::new("blast-task")
                .long("blast-task")
                .default_value("blastn")
                .value_parser(["megablast", "dc-megablast", "blastn", "blastn-short"]),
        )
        .arg(
            Arg::new("evalue")
                .long("evalue")
                .default_value("1e-5")
                .value_parser(value_parser!(f64)),
        )
        .arg(Arg::new("aligner").long("aligner").default_value("mafft"))
        .arg(
            Arg::new("diversity-filter")
                .long("diversity-filter")
                .help("Skip selected loci starting within 100 bp of an already selected locus.")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("force")
                .long("force")
                .action(ArgAction::SetTrue)
                .conflicts_with("resume"),
        )
        .arg(Arg::new("resume").long("resume").action(ArgAction::SetTrue))
}

pub fn parse_recover_copies_args() -> Result<crate::CopyRecoveryConfig> {
    // argv[0] is the binary and argv[1] is the `recover-copies` subcommand name
    recover_copies_from(std::env::args_os().skip(1))
}

fn recover_copies_from<I, T>(args: I) -> Result<crate::CopyRecoveryConfig>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let m = recover_copies_command()
        .try_get_matches_from(args)
        .unwrap_or_else(|e| e.exit());
    let path = |n: &str| m.get_one::<PathBuf>(n).cloned().expect("required by clap");
    let num = |n: &str| *m.get_one::<usize>(n).expect("defaulted by clap");
    let config = crate::CopyRecoveryConfig {
        genome: path("genome"),
        queries: path("queries"),
        output_dir: path("out"),
        flank_length: num("flank-length"),
        hits_per_family: num("hits-per-family"),
        max_fragment_gap: num("max-fragment-gap"),
        orientation_threshold: *m
            .get_one::<f64>("orientation-threshold")
            .expect("defaulted"),
        threads: num("threads"),
        search_backend: m
            .get_one::<String>("search-backend")
            .expect("defaulted")
            .parse()
            .map_err(|e: String| Error::new(ErrorKind::GenericCli(e)))?,
        blast_task: m
            .get_one::<String>("blast-task")
            .expect("defaulted")
            .clone(),
        evalue: *m.get_one::<f64>("evalue").expect("defaulted"),
        aligner: m
            .get_one::<String>("aligner")
            .expect("defaulted")
            .parse()
            .map_err(|e: String| Error::new(ErrorKind::GenericCli(e)))?,
        force: m.get_flag("force"),
        resume: m.get_flag("resume"),
        diversity_filter: m.get_flag("diversity-filter"),
        reputils: None,
    };
    for (label, p) in [("genome", &config.genome), ("queries", &config.queries)] {
        if !p.is_file() {
            return Err(Error::new(ErrorKind::GenericCli(format!(
                "--{label} is not a regular file: {}",
                p.display()
            ))));
        }
    }
    config
        .validate()
        .map_err(|e| Error::new(ErrorKind::GenericCli(e)))?;
    Ok(config)
}

fn classify_command() -> clap::Command {
    use clap::{Arg, Command};
    Command::new("classify")
        .about("Classify RepeatModeler Unknown families and flag candidate novel TEs.")
        .arg(
            Arg::new("library")
                .long("library")
                .required(true)
                .help("RepeatModeler/RepeatMasker library, e.g. consensi.fa.classified")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("refdb")
                .long("refdb")
                .help("Built `rep refdb` database: uses its protein and nucleotide libraries")
                .value_parser(value_parser!(PathBuf))
                .conflicts_with_all(["proteins", "nucleotides"]),
        )
        .arg(
            Arg::new("proteins")
                .long("proteins")
                .help("TE protein FASTA (optionally .gz), e.g. RepBase/RepeatPeps proteins")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("nucleotides")
                .long("nucleotides")
                .help("TE nucleotide library with `name#Class/Family` IDs (optionally .gz)")
                .value_parser(value_parser!(PathBuf)),
        )
        .group(
            clap::ArgGroup::new("reference")
                .args(["refdb", "proteins", "nucleotides"])
                .required(true)
                .multiple(true),
        )
        .arg(
            Arg::new("out")
                .long("out")
                .required(true)
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("all")
                .long("all")
                .help("Classify every family, not only Unknowns, and report label agreement")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .default_value("1")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("evalue")
                .long("evalue")
                .default_value("1e-5")
                .value_parser(value_parser!(f64)),
        )
        .arg(
            Arg::new("backend")
                .long("backend")
                .default_value("blastx")
                .value_parser(["blastx", "diamond"]),
        )
        .arg(
            Arg::new("rm-out")
                .long("rm-out")
                .help("RepeatMasker .out of the genome with this library: copy number, age and genome bp")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("genome")
                .long("genome")
                .requires("rm-out")
                .help("Genome FASTA (plain or .gz) for TSD and flank analysis with --rm-out")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("nucl-task")
                .long("nucl-task")
                .default_value("blastn")
                .value_parser(["blastn", "dc-megablast", "megablast"]),
        )
        .arg(
            Arg::new("nucl-evalue")
                .long("nucl-evalue")
                .default_value("1e-5")
                .value_parser(value_parser!(f64)),
        )
        .arg(
            Arg::new("copies")
                .long("copies")
                .help("`rep recover-copies` output directory: adds copy number, TSD and flank evidence")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("model")
                .long("model")
                .help("Learned classifier run beside the rule cascade (default with --refdb: <db>/build/classifier.json)")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("relatives")
                .long("relatives")
                .help("TE-called families of related genomes (`rep relatives`): Unknowns nothing else \
                       explains take the label of the best match")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("no-library-relatives")
                .long("no-library-relatives")
                .help("Do not search uncalled Unknowns against this library's own TE-called families")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("te-hmms")
                .long("te-hmms")
                .help("Pressed superfamily profile HMMs (default with --refdb: <db>/hmms/te_hmms.hmm)")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("pfam")
                .long("pfam")
                .help("Full Pfam-A.hmm; TE domain models are extracted for ORF domain architecture")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("host-domains")
                .long("host-domains")
                .help("With --pfam: scan ORFs lacking TE evidence against all of Pfam to flag multicopy host genes")
                .requires("pfam")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("min-orf-aa")
                .long("min-orf-aa")
                .default_value("100")
                .value_parser(value_parser!(usize)),
        )
}

pub fn parse_classify_args() -> Result<crate::classify::ClassifyConfig> {
    // argv[0] is the binary and argv[1] is the `classify` subcommand name
    classify_from(std::env::args_os().skip(1))
}

fn classify_from<I, T>(args: I) -> Result<crate::classify::ClassifyConfig>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let m = classify_command()
        .try_get_matches_from(args)
        .unwrap_or_else(|e| e.exit());
    let path = |n: &str| m.get_one::<PathBuf>(n).cloned().expect("required by clap");
    let mut config = if let Some(db) = m.get_one::<PathBuf>("refdb") {
        crate::classify::ClassifyConfig::from_refdb(
            path("library"),
            &crate::refdb::RefDb::new(db.clone()),
            path("out"),
        )
    } else {
        let mut c =
            crate::classify::ClassifyConfig::new(path("library"), PathBuf::new(), path("out"));
        c.proteins = m.get_one::<PathBuf>("proteins").cloned();
        c.nucleotides = m.get_one::<PathBuf>("nucleotides").cloned();
        c
    };
    config.all_families = m.get_flag("all");
    config.pfam = m.get_one::<PathBuf>("pfam").cloned();
    config.host_domains = m.get_flag("host-domains");
    if let Some(hmms) = m.get_one::<PathBuf>("te-hmms") {
        config.te_hmms = Some(hmms.clone());
    }
    if let Some(model) = m.get_one::<PathBuf>("model") {
        config.model = Some(model.clone());
    }
    config.relatives = m.get_one::<PathBuf>("relatives").cloned();
    config.library_relatives = !m.get_flag("no-library-relatives");
    config.copies_dir = m.get_one::<PathBuf>("copies").cloned();
    config.rm_out = m.get_one::<PathBuf>("rm-out").cloned();
    config.genome = m.get_one::<PathBuf>("genome").cloned();
    config.nucl_task = m.get_one::<String>("nucl-task").expect("defaulted").clone();
    config.nucl_evalue = *m.get_one::<f64>("nucl-evalue").expect("defaulted");
    config.threads = *m.get_one::<usize>("threads").expect("defaulted");
    config.evalue = *m.get_one::<f64>("evalue").expect("defaulted");
    config.backend = m
        .get_one::<String>("backend")
        .expect("defaulted")
        .parse()
        .map_err(|e: String| Error::new(ErrorKind::GenericCli(e)))?;
    config.structure.min_orf_aa = *m.get_one::<usize>("min-orf-aa").expect("defaulted");
    config
        .validate()
        .map_err(|e| Error::new(ErrorKind::GenericCli(e)))?;
    Ok(config)
}

fn refdb_command() -> clap::Command {
    use clap::{Arg, Command};
    let db = || {
        Arg::new("db")
            .long("db")
            .required(true)
            .help("Reference database directory")
            .value_parser(value_parser!(PathBuf))
    };
    Command::new("refdb")
        .about("Manage the merged Repbase + Dfam + local TE reference database.")
        .subcommand_required(true)
        .subcommand(
            Command::new("fetch")
                .about("Download and register the latest (or a given) public releases")
                .arg(db())
                .arg(
                    Arg::new("repbase-release")
                        .long("repbase-release")
                        .help("e.g. RepBase31.09"),
                )
                .arg(
                    Arg::new("no-repbase")
                        .long("no-repbase")
                        .action(ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("no-dfam")
                        .long("no-dfam")
                        .action(ArgAction::SetTrue),
                ),
        )
        .subcommand(
            Command::new("register")
                .about("Register an already-downloaded source release (offline use)")
                .arg(db())
                .arg(
                    Arg::new("source")
                        .long("source")
                        .required(true)
                        .value_parser(["repbase", "dfam-curated", "dfam-peptides"]),
                )
                .arg(Arg::new("version").long("version").required(true))
                .arg(
                    Arg::new("file")
                        .long("file")
                        .required(true)
                        .value_parser(value_parser!(PathBuf)),
                ),
        )
        .subcommand(
            Command::new("add")
                .about("Add our own families to the local layer")
                .arg(db())
                .arg(
                    Arg::new("fasta")
                        .long("fasta")
                        .value_parser(value_parser!(PathBuf)),
                )
                .arg(
                    Arg::new("from-classify")
                        .long("from-classify")
                        .help("A `rep classify` output directory")
                        .value_parser(value_parser!(PathBuf))
                        .conflicts_with("fasta"),
                )
                .arg(
                    Arg::new("dispositions")
                        .long("dispositions")
                        .default_value("novel_candidate")
                        .help("Comma-separated dispositions taken from --from-classify"),
                )
                .arg(
                    Arg::new("kind")
                        .long("kind")
                        .default_value("nucleotide")
                        .value_parser(["nucleotide", "protein"]),
                )
                .arg(
                    Arg::new("origin")
                        .long("origin")
                        .required(true)
                        .help("Genome or project tag, prefixed to IDs (e.g. Ec32)"),
                )
                .arg(
                    Arg::new("status")
                        .long("status")
                        .default_value("candidate")
                        .value_parser(["accepted", "candidate", "rejected"]),
                )
                .arg(Arg::new("notes").long("notes").default_value("")),
        )
        .subcommand(
            Command::new("build")
                .about("Merge sources into build/")
                .arg(db()),
        )
        .subcommand(
            Command::new("hmms")
                .about("Build superfamily profile HMMs from the built protein library (hmms/te_hmms.hmm)")
                .arg(db())
                .arg(
                    Arg::new("threads")
                        .long("threads")
                        .default_value("1")
                        .value_parser(value_parser!(usize)),
                ),
        )
        .subcommand(
            Command::new("mosaics")
                .about("Find references carrying an embedded element of another class (build/mosaics.tsv)")
                .arg(db())
                .arg(
                    Arg::new("threads")
                        .long("threads")
                        .default_value("1")
                        .value_parser(value_parser!(usize)),
                ),
        )
        .subcommand(
            Command::new("status")
                .about("Show sources, local layer and build")
                .arg(db()),
        )
}

/// Parse and run `rep refdb ...`.
pub fn run_refdb() -> Result<()> {
    use crate::refdb;
    let m = refdb_command()
        .try_get_matches_from(std::env::args_os().skip(1))
        .unwrap_or_else(|e| e.exit());
    let (name, sub) = m.subcommand().expect("subcommand_required");
    let db = refdb::RefDb::new(sub.get_one::<PathBuf>("db").expect("required").clone());
    let runner = crate::RealCommandRunner;
    let text = |n: &str| sub.get_one::<String>(n).cloned().unwrap_or_default();
    match name {
        "fetch" => {
            let options = refdb::FetchOptions {
                repbase: !sub.get_flag("no-repbase"),
                repbase_release: sub.get_one::<String>("repbase-release").cloned(),
                dfam: !sub.get_flag("no-dfam"),
            };
            let added = refdb::fetch(&db, &options, &runner)?;
            eprintln!(
                "Registered {} new release(s); run `rep refdb build` to merge",
                added.len()
            );
        }
        "register" => {
            let file = sub.get_one::<PathBuf>("file").expect("required");
            let release = refdb::register_source(
                &db,
                &text("source"),
                &text("version"),
                file,
                None,
                &runner,
            )?;
            eprintln!("Registered {} {}", release.name, release.version);
        }
        "add" => {
            let status = text("status");
            let records = if let Some(dir) = sub.get_one::<PathBuf>("from-classify") {
                let dispositions = text("dispositions")
                    .split(',')
                    .map(str::to_owned)
                    .collect::<Vec<_>>();
                refdb::local_records_from_classify(dir, &dispositions, &text("origin"), &status)?
            } else if let Some(fasta) = sub.get_one::<PathBuf>("fasta") {
                let kind = text("kind")
                    .parse()
                    .map_err(|e: String| Error::new(ErrorKind::GenericCli(e)))?;
                refdb::local_records_from_fasta(
                    fasta,
                    kind,
                    &text("origin"),
                    &status,
                    &text("notes"),
                )?
            } else {
                return Err(Error::new(ErrorKind::GenericCli(
                    "refdb add needs --fasta or --from-classify".into(),
                )));
            };
            let summary = db.add_local(records)?;
            eprintln!(
                "Added {} local record(s), {} already present",
                summary.added, summary.unchanged
            );
        }
        "mosaics" => {
            let options = refdb::MosaicOptions {
                threads: *sub.get_one::<usize>("threads").expect("defaulted"),
                ..Default::default()
            };
            let info = refdb::find_mosaics(&db, &options, &runner)?;
            eprintln!(
                "{} of {} references carry {} embedded segment(s): {:?}",
                info.mosaic_references, info.references, info.segments, info.by_pair
            );
        }
        "hmms" => {
            let options = refdb::HmmBuildOptions {
                threads: *sub.get_one::<usize>("threads").expect("defaulted"),
                ..Default::default()
            };
            let info = refdb::build_te_hmms(&db, &options, &runner)?;
            eprintln!(
                "Built {} profile HMMs for {} labels from {} proteins ({} clusters, {} failed)",
                info.models, info.labels, info.proteins, info.clusters, info.failed
            );
        }
        "build" => {
            let info = refdb::build(&db)?;
            eprintln!(
                "Built {} nucleotide and {} protein records ({} duplicates, {} renamed, {} unlabelled)",
                info.nucleotide_records,
                info.protein_records,
                info.duplicates,
                info.renamed,
                info.unlabelled
            );
        }
        _ => print!("{}", refdb::status_report(&db)?),
    }
    Ok(())
}

fn discover_command() -> clap::Command {
    use clap::{Arg, Command};
    let num = |name: &'static str, default: &'static str| {
        Arg::new(name)
            .long(name)
            .default_value(default)
            .value_parser(value_parser!(f64))
    };
    Command::new("discover")
        .about("Cluster families across many `rep classify --all` runs; propagate labels and rank novel groups.")
        .arg(
            Arg::new("results")
                .long("results")
                .required(true)
                .help("Directory searched for classify outputs (classification.tsv + queries.fa)")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(Arg::new("out").long("out").required(true).value_parser(value_parser!(PathBuf)))
        .arg(
            Arg::new("taxa")
                .long("taxa")
                .help("TSV with scientific_name and class columns")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(Arg::new("threads").long("threads").default_value("1").value_parser(value_parser!(usize)))
        .arg(num("family-identity", "0.8"))
        .arg(num("family-coverage", "0.8"))
        .arg(num("protein-identity", "0.3"))
        .arg(num("protein-coverage", "0.5"))
        .arg(num("min-purity", "0.75"))
}

/// Parse and run `rep relatives ...`.
pub fn run_relatives() -> Result<()> {
    use clap::{Arg, Command};
    let m = Command::new("relatives")
        .about("Pool the families `rep classify` runs called as TEs, one FASTA per taxon, for `classify --relatives`.")
        .arg(
            Arg::new("results")
                .long("results")
                .required(true)
                .help("Directory searched for classify outputs (classification.tsv + queries.fa)")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(Arg::new("out").long("out").required(true).value_parser(value_parser!(PathBuf)))
        .arg(
            Arg::new("taxa")
                .long("taxa")
                .help("TSV with a scientific_name (or species) column and the --rank column; without it, one file")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(Arg::new("rank").long("rank").default_value("order").help("Taxonomy column that groups genomes"))
        .try_get_matches_from(std::env::args_os().skip(1))
        .unwrap_or_else(|e| e.exit());
    let path = |n: &str| m.get_one::<PathBuf>(n).cloned().expect("required by clap");
    let s = crate::discover::write_relatives(&crate::discover::RelativesConfig {
        results: path("results"),
        taxa: m.get_one::<PathBuf>("taxa").cloned(),
        rank: m.get_one::<String>("rank").expect("defaulted").clone(),
        output_dir: path("out"),
    })?;
    eprintln!(
        "{} reference families from {} samples in {} taxa",
        s.per_taxon.values().sum::<usize>(),
        s.samples,
        s.per_taxon.len()
    );
    Ok(())
}

/// Parse and run `rep discover ...`.
pub fn run_discover() -> Result<()> {
    let m = discover_command()
        .try_get_matches_from(std::env::args_os().skip(1))
        .unwrap_or_else(|e| e.exit());
    let path = |n: &str| m.get_one::<PathBuf>(n).cloned().expect("required by clap");
    let num = |n: &str| *m.get_one::<f64>(n).expect("defaulted");
    let mut c = crate::discover::DiscoverConfig::new(path("results"), path("out"));
    c.taxa = m.get_one::<PathBuf>("taxa").cloned();
    c.threads = *m.get_one::<usize>("threads").expect("defaulted");
    c.family_identity = num("family-identity");
    c.family_coverage = num("family-coverage");
    c.protein_identity = num("protein-identity");
    c.protein_coverage = num("protein-coverage");
    c.min_purity = num("min-purity");
    let s = crate::discover::run_discover(&c, &crate::RealCommandRunner)?;
    eprintln!(
        "{} families from {} samples: {} family clusters, {} groups; propagated {} (family) + {} (group); \
         leave-one-out accuracy {:.3} (family, n={}) / {:.3} (group, n={}); {} novel groups",
        s.members,
        s.samples,
        s.family_clusters,
        s.groups,
        s.propagated_from_family_cluster,
        s.propagated_from_group,
        s.loo_family_cluster.accuracy,
        s.loo_family_cluster.predicted,
        s.loo_group.accuracy,
        s.loo_group.predicted,
        s.novel_groups
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recover_copies_parses_typed_options() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("x.fa");
        std::fs::write(&fa, ">a\nACGT\n").unwrap();
        let fa = fa.to_str().unwrap();
        let c = recover_copies_from([
            "recover-copies",
            "--queries",
            fa,
            "--genome",
            fa,
            "--out",
            "o",
            "--evalue",
            "1e-3",
            "--blast-task",
            "dc-megablast",
            "--diversity-filter",
        ])
        .unwrap();
        assert_eq!(c.blast_task, "dc-megablast");
        assert_eq!(c.evalue, 1e-3);
        assert!(c.diversity_filter);
        assert_eq!(c.hits_per_family, 20);
        let d = recover_copies_from([
            "recover-copies",
            "--queries",
            fa,
            "--genome",
            fa,
            "--out",
            "o",
        ])
        .unwrap();
        assert_eq!(d.blast_task, crate::copy_recovery::default_blast_task());
        assert_eq!(d.evalue, crate::copy_recovery::default_evalue());
        assert_eq!((d.blast_task.as_str(), d.evalue), ("blastn", 1e-5));
    }

    #[test]
    fn recover_copies_rejects_missing_inputs() {
        let err = recover_copies_from([
            "recover-copies",
            "--queries",
            "/nonexistent.fa",
            "--genome",
            "/nonexistent.fa",
            "--out",
            "o",
        ])
        .unwrap_err();
        assert!(err.to_string().contains("not a regular file"));
    }

    #[test]
    fn classify_parses_options() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("x.fa");
        std::fs::write(&fa, ">a\nACGT\n").unwrap();
        let fa = fa.to_str().unwrap();
        let c = classify_from([
            "classify",
            "--library",
            fa,
            "--proteins",
            fa,
            "--out",
            "o",
            "--all",
            "--backend",
            "diamond",
            "--min-orf-aa",
            "80",
        ])
        .unwrap();
        assert!(c.all_families);
        assert_eq!(c.backend, crate::classify::HomologyBackend::Diamond);
        assert_eq!(c.structure.min_orf_aa, 80);
    }
}
