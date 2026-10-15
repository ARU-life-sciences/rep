use super::*;
use std::os::unix::process::ExitStatusExt;
use std::process::{ExitStatus, Output};

/// Fakes BLAST: writes fixed hits for the reference searches.
struct MockBlast {
    nucl: String,
    prot: String,
}

impl CommandRunner for MockBlast {
    fn run(&self, cmd: &mut Command) -> Result<Output> {
        let args = cmd
            .get_args()
            .map(|a| a.to_string_lossy().into_owned())
            .collect::<Vec<_>>();
        if let Some(i) = args.iter().position(|a| a == "-out") {
            let out = &args[i + 1];
            if out.ends_with("nucleotide_hits.0.tsv") {
                fs::write(out, &self.nucl)?;
            } else if out.ends_with("protein_hits.tsv") {
                fs::write(out, &self.prot)?;
            }
        }
        Ok(Output {
            status: ExitStatus::from_raw(0),
            stdout: vec![],
            stderr: vec![],
        })
    }
}

fn random(n: usize, seed: u64) -> String {
    let mut x = seed;
    (0..n)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ['A', 'C', 'G', 'T'][(x >> 62) as usize]
        })
        .collect()
}

#[test]
fn hits_resolve_chain_and_mask_with_tandem_runs_in_the_gaps() {
    let dir = tempfile::tempdir().unwrap();
    let db = dir.path().join("db");
    fs::create_dir_all(db.join("build")).unwrap();
    fs::write(db.join("build/te_nucl.fa"), ">x#DNA\nACGT\n").unwrap();
    fs::write(db.join("build/te_prot.fa"), ">p#DNA\nMK\n").unwrap();
    let ctg = format!("{}{}{}", random(2000, 1), "CA".repeat(40), random(920, 2));
    let input = dir.path().join("in.fa");
    fs::write(
        &input,
        format!(">ctg1 a contig\n{ctg}\n>ctg2\n{}\n", random(500, 3)),
    )
    .unwrap();
    let row = |q: &str, s: &str, slen, qs, qe, ss, se, bits| {
        format!("{q}\t{s}\t3000\t{slen}\t85.0\t500\t{qs}\t{qe}\t{ss}\t{se}\t1e-50\t{bits}\n")
    };
    let runner = MockBlast {
        nucl: [
            row("q0_0", "Gypsy-1#LTR/Gypsy", 5000, 101, 600, 1, 500, 500.0),
            // overlaps the Gypsy hit, on the minus strand: keeps 601..900
            row("q0_0", "L1-2#LINE/L1", 6000, 400, 900, 900, 400, 300.0),
        ]
        .concat(),
        // the first unclaimed stretch long enough for proteins starts at 900
        prot: row("q0_900", "RTE-1_p#LINE/RTE", 900, 101, 400, 10, 109, 120.0),
    };
    let mut c = AnnotateConfig::new(input, dir.path().join("out"));
    c.refdb = Some(db);
    let s = run_annotate(&c, &runner).unwrap();
    let table = fs::read_to_string(dir.path().join("out").join(TABLE)).unwrap();
    let rows = table
        .lines()
        .skip(1)
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let cols = |r: &Vec<&str>| {
        (
            r[1].to_owned(),
            r[2].to_owned(),
            r[4].to_owned(),
            r[12].to_owned(),
        )
    };
    assert_eq!(
        cols(&rows[0]),
        (
            "101".into(),
            "600".into(),
            "LTR/Gypsy".into(),
            "nucleotide".into()
        )
    );
    // the trimmed minus-strand hit keeps the matching end of its reference
    assert_eq!(cols(&rows[1]).0, "601");
    assert_eq!((rows[1][3], rows[1][6], rows[1][7]), ("-", "400", "699"));
    assert_eq!(
        cols(&rows[2]),
        (
            "1001".into(),
            "1300".into(),
            "LINE/RTE".into(),
            "protein".into()
        )
    );
    assert_eq!((rows[3][4], rows[3][12]), ("Simple_repeat", "tandem"));
    assert!(
        rows[3][5] == "(CA)n" || rows[3][5] == "(AC)n",
        "{}",
        rows[3][5]
    );
    assert_eq!(rows.len(), 4, "{table}");
    assert_eq!(s.elements, 4);
    assert_eq!((s.by_class["LTR"], s.by_class["LINE"]), (500, 600));
    assert!(s.hint.is_none());
    let masked = fs::read_to_string(dir.path().join("out").join(MASKED)).unwrap();
    let seq1 = masked.lines().skip(1).take(50).collect::<String>();
    assert_eq!(seq1.len(), 3000);
    assert!(seq1[100..600].bytes().all(|b| b.is_ascii_lowercase()));
    assert!(seq1[..100].bytes().all(|b| b.is_ascii_uppercase()));
    let gff = fs::read_to_string(dir.path().join("out").join(GFF)).unwrap();
    assert!(gff.contains(
        "ctg1\trep\ttransposable_element\t101\t600\t500.0\t+\t.\tID=frag1;Name=Gypsy-1#LTR/Gypsy;class=LTR/Gypsy;element=E1;evidence=nucleotide;Target=Gypsy-1#LTR/Gypsy 1 500"
    ));
    assert!(gff.contains("\tmicrosatellite\t"));
}

#[test]
fn one_short_sequence_covered_by_one_element_gets_a_hint() {
    let dir = tempfile::tempdir().unwrap();
    let lib = dir.path().join("lib.fa");
    fs::write(&lib, ">fam1#DNA/hAT\nACGT\n").unwrap();
    let input = dir.path().join("in.fa");
    fs::write(&input, format!(">element\n{}\n", random(1000, 9))).unwrap();
    let runner = MockBlast {
        nucl: "q0_0\tfam1#DNA/hAT\t1000\t1000\t99.0\t950\t21\t970\t1\t950\t0\t1700\n".into(),
        prot: String::new(),
    };
    let mut c = AnnotateConfig::new(input, dir.path().join("out"));
    c.libraries = vec![lib];
    let s = run_annotate(&c, &runner).unwrap();
    assert!(s.hint.unwrap().contains("rep classify"));
    assert_eq!(s.by_class["DNA"], 950);
    let table = fs::read_to_string(dir.path().join("out").join(TABLE)).unwrap();
    assert!(table.contains("\tlib.fa\t"));
}

#[test]
fn windows_overlap_and_cover_the_stretch() {
    assert_eq!(windows(100, 250, 100, 20), vec![(100, 200), (180, 250)]);
    assert_eq!(windows(0, 50, 100, 20), vec![(0, 50)]);
    assert_eq!(window_of("q3_1800").unwrap(), (3, 1800));
    assert!(window_of("x").is_err());
}

#[test]
fn classified_references_claim_bases_before_unknown_ones() {
    let hit = |reference: &str, class: &str, bitscore: f64| Fragment {
        start: 0,
        end: 1000,
        reference: reference.into(),
        class: class.into(),
        bitscore,
        ..Default::default()
    };
    let order = by_score(vec![
        hit("DM412#Unknown", "Unknown", 2000.0),
        hit("DM412_I#LTR/Gypsy", "LTR/Gypsy", 1500.0),
        hit("DM412_LTR#LTR/Gypsy", "LTR/Gypsy", 1600.0),
    ]);
    assert_eq!(
        order
            .iter()
            .map(|f| f.reference.as_str())
            .collect::<Vec<_>>(),
        vec!["DM412_LTR#LTR/Gypsy", "DM412_I#LTR/Gypsy", "DM412#Unknown"]
    );
}
