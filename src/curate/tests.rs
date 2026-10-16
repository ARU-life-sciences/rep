use super::*;

fn random(n: usize, seed: u64) -> Vec<u8> {
    let mut x = seed;
    (0..n)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(x >> 62) as usize]
        })
        .collect()
}

/// `seq` with every `every`-th base (from `offset`) changed.
fn mutate(seq: &[u8], every: usize, offset: usize) -> Vec<u8> {
    seq.iter()
        .enumerate()
        .map(|(i, &b)| {
            if i % every == offset {
                match b {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                }
            } else {
                b
            }
        })
        .collect()
}

/// Writes a recover-copies family directory whose copies are already aligned
/// (no indels): each copy is `flank + element + flank`, and its BLAST core is
/// `element[core.0..core.1]`.
fn family(dir: &Path, query: &[u8], copies: &[Vec<u8>], flank: usize, core: (usize, usize)) {
    fs::create_dir_all(dir).unwrap();
    fs::write(
        dir.join("query.fa"),
        format!(">fam1#Unknown\n{}\n", String::from_utf8_lossy(query)),
    )
    .unwrap();
    let mut aligned = String::new();
    let mut tsv = String::from("family_id\tlocus_id\tsubject_id\tgenomic_start\tgenomic_end\tstrand\trequested_extended_start\trequested_extended_end\tactual_extended_start\tactual_extended_end\tleft_flank_actual\tright_flank_actual\treverse_complemented\tselected\textraction_status\talignment_status\texclusion_reason\n");
    let mut loci = String::from("family_id\tlocus_id\tsubject_id\tstart\tend\textended_start\textended_end\tstrand\tn_fragments\tbest_bitscore\tcombined_query_coverage\tpercent_identity\trepresentative_score\tselected\texclusion_reason\n");
    for (i, c) in copies.iter().enumerate() {
        let id = format!("locus-{i:06}");
        aligned += &format!(">locus_id={id} family_id=fam1#Unknown\n{}\n", String::from_utf8_lossy(c));
        let left = flank + core.0;
        let right = c.len() - flank - core.1;
        let es = 10_001 + 100_000 * i;
        let ee = es + c.len() - 1;
        tsv += &format!(
            "fam1#Unknown\t{id}\tchr1\t{}\t{}\t+\t{es}\t{ee}\t{es}\t{ee}\t{left}\t{right}\tfalse\ttrue\tsuccess\tpending\t\n",
            es + left,
            ee - right
        );
        loci += &format!("fam1#Unknown\t{id}\tchr1\t0\t0\t0\t0\t+\t1\t500\t1.0\t95.0\t1\ttrue\t\n");
    }
    fs::write(dir.join("copies.aligned.fa"), &aligned).unwrap();
    fs::write(dir.join("copies.extended.fa"), aligned.replace('-', "")).unwrap();
    fs::write(dir.join("copies.tsv"), tsv).unwrap();
    fs::write(dir.join("loci.tsv"), loci).unwrap();
}

fn element_copies(element: &[u8], n: usize, flank: usize) -> Vec<Vec<u8>> {
    (0..n)
        .map(|i| {
            [
                random(flank, 100 + i as u64),
                mutate(element, 17, i),
                random(flank, 200 + i as u64),
            ]
            .concat()
        })
        .collect()
}

#[test]
fn a_fragment_query_is_extended_to_the_element_ends() {
    let dir = tempfile::tempdir().unwrap();
    let element = random(600, 1);
    let rec = dir.path().join("rec");
    family(
        &rec.join("families/fam1_Unknown"),
        &element[150..450],
        &element_copies(&element, 8, 300),
        300,
        (150, 450),
    );
    let c = CurateConfig::new(rec, dir.path().join("out"));
    let s = run_curate(&c).unwrap();
    assert_eq!((s.families, s.refined), (1, 1));
    let refined = read_fasta(&dir.path().join("out").join(REFINED)).unwrap();
    assert_eq!(refined[0].0, "fam1#Unknown");
    // flank columns can agree by chance for a few bases past the end
    let r = &refined[0].1;
    let at = r.windows(element.len()).position(|w| w == element).expect("whole element");
    assert!(at <= 5 && r.len() - element.len() <= 10, "{at} {}", r.len());
    let table = fs::read_to_string(dir.path().join("out").join(TABLE)).unwrap();
    assert!(table.contains("fam1#Unknown\t300\t"), "{table}");
    assert!(table.trim_end().ends_with("\t8\trefined"), "{table}");
    // copy boundaries now sit at the element ends, for TSDs and full-length counts
    let tsv = fs::read_to_string(dir.path().join("out/families/fam1_Unknown/copies.tsv")).unwrap();
    let first = tsv.lines().nth(1).unwrap().split('\t').collect::<Vec<_>>();
    let (left, right): (usize, usize) = (first[10].parse().unwrap(), first[11].parse().unwrap());
    assert!(left.abs_diff(300) <= 5 && right.abs_diff(300) <= 10, "{left} {right}");
    assert_eq!(first[3], (10_001 + left).to_string());
    let loci = fs::read_to_string(dir.path().join("out/families/fam1_Unknown/loci.tsv")).unwrap();
    let cov: f64 = loci.lines().nth(1).unwrap().split('\t').nth(10).unwrap().parse().unwrap();
    assert!(cov >= 0.95, "{loci}");
}

#[test]
fn copies_ending_together_stop_the_extension() {
    // half the copies continue into shared sequence (LTRs of full elements),
    // half into unrelated flanks (solo LTRs): the boundary is the LTR end
    let ltr = random(300, 2);
    let internal = random(300, 3);
    let copies = (0..8)
        .map(|i| {
            let right = if i % 2 == 0 {
                mutate(&internal, 13, i)
            } else {
                random(300, 300 + i as u64)
            };
            [random(300, 100 + i as u64), mutate(&ltr, 17, i), right].concat()
        })
        .collect::<Vec<_>>();
    let rows = copies
        .iter()
        .enumerate()
        .map(|(i, s)| Row {
            locus_id: i.to_string(),
            seq: s.clone(),
            first: 0,
            last: s.len() - 1,
            core: (300, 599),
        })
        .collect::<Vec<_>>();
    let r = refine(&rows, 300, &CurateParams::default()).unwrap();
    assert_eq!(r.seq, ltr);
}

#[test]
fn minority_insertions_are_dropped_not_written_as_n() {
    let element = random(400, 4);
    let mut copies = element_copies(&element, 6, 100)
        .into_iter()
        .map(|c| {
            let mut c = c;
            c.insert(300, b'-');
            c
        })
        .collect::<Vec<_>>();
    copies[0][300] = b'G';
    let rows = copies
        .iter()
        .enumerate()
        .map(|(i, s)| Row {
            locus_id: i.to_string(),
            seq: s.clone(),
            first: 0,
            last: s.len() - 1,
            core: (100, 500),
        })
        .collect::<Vec<_>>();
    let r = refine(&rows, 400, &CurateParams::default()).unwrap();
    assert_eq!(r.seq, element);
}

#[test]
fn too_few_copies_or_a_short_result_keep_the_query() {
    let element = random(600, 5);
    let rows = |n: usize, core: (usize, usize)| {
        element_copies(&element, n, 300)
            .into_iter()
            .enumerate()
            .map(|(i, s)| Row {
                locus_id: i.to_string(),
                first: 0,
                last: s.len() - 1,
                seq: s,
                core,
            })
            .collect::<Vec<_>>()
    };
    let p = CurateParams::default();
    assert_eq!(refine(&rows(2, (300, 899)), 600, &p), Err("too_few_copies"));
    // a query twice the length the copies support
    assert_eq!(refine(&rows(5, (300, 899)), 1500, &p), Err("shorter_than_query"));

    let dir = tempfile::tempdir().unwrap();
    let rec = dir.path().join("rec");
    family(
        &rec.join("families/fam1_Unknown"),
        &element,
        &element_copies(&element, 2, 300),
        300,
        (0, 600),
    );
    let s = run_curate(&CurateConfig::new(rec, dir.path().join("out"))).unwrap();
    assert_eq!((s.refined, s.kept["too_few_copies"]), (0, 1));
    let refined = read_fasta(&dir.path().join("out").join(REFINED)).unwrap();
    assert_eq!(refined[0].1, element);
    assert!(dir.path().join("out/families/fam1_Unknown/copies.tsv").is_file());
}

#[test]
fn extended_ends_snap_to_the_elements_own_inverted_repeats() {
    // a TIR element whose copies all sit beside the same 120 bp of another
    // repeat, so they agree past its 3' end; its TIRs mark the real end
    let tir = random(30, 6);
    let body = random(400, 7);
    let element = [tir.clone(), body, bio::alphabets::dna::revcomp(&tir)].concat();
    let neighbour = random(120, 8);
    let copies = (0..8)
        .map(|i| {
            [
                random(300, 100 + i as u64),
                mutate(&element, 23, i),
                mutate(&neighbour, 19, i),
                random(300, 200 + i as u64),
            ]
            .concat()
        })
        .collect::<Vec<_>>();
    let rows = copies
        .iter()
        .enumerate()
        .map(|(i, s)| Row {
            locus_id: i.to_string(),
            seq: s.clone(),
            first: 0,
            last: s.len() - 1,
            core: (400, 659),
        })
        .collect::<Vec<_>>();
    let r = refine(&rows, 260, &CurateParams::default()).unwrap();
    assert_eq!(r.seq, element);
    assert_eq!((r.added_5, r.added_3), (100, 100));
}

#[test]
fn a_conserved_tsd_taken_into_the_extension_is_trimmed_off() {
    // a Tc1-like element: every copy sits in TA..TA, so the extension takes
    // the TSD (and the copies' agreement on it) into the consensus
    let element = [b"CAGT".to_vec(), random(600, 9), b"ACTG".to_vec()].concat();
    let copies = (0..12)
        .map(|i| {
            [
                random(300, 100 + i as u64),
                b"TA".to_vec(),
                mutate(&element, 23, i),
                b"TA".to_vec(),
                random(300, 200 + i as u64),
            ]
            .concat()
        })
        .collect::<Vec<_>>();
    let rows = copies
        .iter()
        .enumerate()
        .map(|(i, s)| Row {
            locus_id: i.to_string(),
            seq: s.clone(),
            first: 0,
            last: s.len() - 1,
            core: (402, 701),
        })
        .collect::<Vec<_>>();
    let r = refine(&rows, 300, &CurateParams::default()).unwrap();
    assert_eq!(r.seq, element);
}

#[test]
fn the_review_package_holds_both_consensi_and_the_end_views() {
    let dir = tempfile::tempdir().unwrap();
    let element = random(600, 1);
    let rec = dir.path().join("rec");
    family(
        &rec.join("families/fam1_Unknown"),
        &element[150..450],
        &element_copies(&element, 8, 300),
        300,
        (150, 450),
    );
    run_curate(&CurateConfig::new(rec, dir.path().join("out"))).unwrap();
    let fam = dir.path().join("out/families/fam1_Unknown");
    let aligned = read_fasta(&fam.join("review.aligned.fa")).unwrap();
    assert_eq!(aligned.len(), 2 + 8);
    assert_eq!(aligned[0].0, "refined_consensus");
    assert_eq!(aligned[1].0, "raw_consensus");
    assert!(aligned.iter().all(|(_, s)| s.len() == aligned[2].1.len()));
    // the raw consensus sits in the columns of the copies it came from
    let ungapped = |s: &[u8]| s.iter().copied().filter(|&b| b != b'-').collect::<Vec<_>>();
    assert_eq!(ungapped(&aligned[1].1), &element[150..450]);
    let at = aligned[1].1.iter().position(|&b| b != b'-').unwrap();
    // (copies carry a planted mutation every 17 bases)
    let same = (at..at + 34).filter(|&c| aligned[2].1[c] == aligned[1].1[c]).count();
    assert!(same >= 32, "{same}");
    let txt = fs::read_to_string(fam.join("review.txt")).unwrap();
    assert!(txt.contains("raw_bases_unplaced\t0"), "{txt}");
    assert!(txt.contains("5' end: alignment columns") && txt.contains("3' end: alignment columns"));
    let review = fs::read_to_string(dir.path().join("out").join(REVIEW)).unwrap();
    let row = review.lines().nth(1).unwrap().split('\t').collect::<Vec<_>>();
    assert_eq!((row[0], row[1], row[12]), ("fam1#Unknown", "refined", "families/fam1_Unknown"));
    assert_eq!(row.len(), 16);
}

#[test]
fn a_degraded_ttaa_tsd_is_trimmed_off_whole() {
    // piggyBac-like: TTAA duplications intact in only half the copies (old
    // insertions), too few for a 2 bp TSD but enough for a 4 bp one
    let element = [b"CCCT".to_vec(), random(600, 9), b"AGGG".to_vec()].concat();
    let copies = (0..12u64)
        .map(|i| {
            let tsd = if i % 2 == 0 { b"TTAA".to_vec() } else { random(4, 300 + i) };
            let right = if i % 2 == 0 { b"TTAA".to_vec() } else { random(4, 400 + i) };
            [random(300, 100 + i), tsd, mutate(&element, 23, i as usize), right, random(300, 200 + i)].concat()
        })
        .collect::<Vec<_>>();
    let rows = copies
        .iter()
        .enumerate()
        .map(|(i, s)| Row {
            locus_id: i.to_string(),
            seq: s.clone(),
            first: 0,
            last: s.len() - 1,
            core: (404, 703),
        })
        .collect::<Vec<_>>();
    let r = refine(&rows, 300, &CurateParams::default()).unwrap();
    let tsd = r.tsd.expect("a TSD trim");
    assert_eq!((tsd.length, tsd.motif.as_str()), (4, "TTAA"));
    assert_eq!(r.seq, element);
}
