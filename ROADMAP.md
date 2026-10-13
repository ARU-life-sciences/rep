# rep roadmap: exhaustive resolution of Unknown TE families and discovery of novel classes

## Goal

1. **Every** RepeatModeler/RepeatMasker family labelled `Unknown` receives an
   evidence-backed *disposition*. No family stays unexplained without a recorded reason.
2. Families that look like real transposable elements but match no known class are
   surfaced and ranked as **candidate novel TE classes**.

A disposition is one of:

| Disposition | Meaning |
| --- | --- |
| `known_te` | Protein, nucleotide or domain-architecture evidence for a known class; relabelled `Class/Superfamily` |
| `chimeric` | Homology to two or more incompatible TE classes in different regions |
| `simple_repeat` / `satellite` | Tandem periodicity explains the family (period recorded) |
| `low_complexity` | Compositionally biased, non-periodic |
| `novel_candidate` | TE hallmarks (terminal repeats, TSDs, coding capacity, unusual domain architecture) without known homology |
| `interspersed_repeat` | A confirmed dispersed multicopy repeat (unique flanks, many full-length copies) with no TE hallmark yet |
| `host_or_artefact` | Matches an RNA gene or other non-TE reference, or keeps such a label from RepeatClassifier |
| `unresolved` | No evidence either way; the missing evidence is recorded |

## Pipeline at a glance

```text
per genome                                                    across genomes
──────────                                                    ──────────────
genome.fa
  │ 1. RepeatModeler            consensus library (RepeatClassifier labels, mostly Unknown)
  │ 2. copy recovery            all genomic copies of each consensus, with flanks
  │    (or an existing RepeatMasker .out + genome: same evidence, much faster)
  │ 3. rep classify             evidence per family → disposition → relabelled library
  │ 4. RepeatMasker             genome annotated with the relabelled library
  ▼                                                           5. rep discover
classification.tsv, masked genome  ─────────────────────────▶  cluster families of many species,
                                                               propagate labels, rank novel groups
                                                               │
             rep refdb (Repbase + Dfam + our own families) ◀──┘ accepted novel families
             build → profile HMMs → mosaics                     (refdb add --from-classify)
```

Evidence gathered in step 3, per family:

| Evidence | Source | Tells us |
| --- | --- | --- |
| Protein homology | BLASTX vs 78k refdb proteins | class of coding families |
| Nucleotide homology | BLASTN vs 146k refdb consensi, mosaic-aware | class of non-coding/diverged families |
| Domain architecture | Pfam TE subset + 1,712 superfamily profile HMMs | class from domain order (e.g. GAG-PR-RT-RH-INT) |
| Host genes | full Pfam on ORFs lacking TE evidence (`--host-domains`) | multicopy host genes, processed pseudogenes |
| Structure | consensus sequence | LTRs, TIRs, poly-A tails, tandem period, termini |
| Copies | genome | copy number, age, TSD length/motif, Helitron ends, flank conservation |

## Toward one classifier (M8, planned)

Today `decide.rs` is a hand-ordered rule cascade: non-TE nucleotide matches first, then protein,
domains, nucleotide, structure, with copy evidence adjusting the result. Confidence bands were
calibrated afterwards (high 98.5%, medium 97.1%, low 84.1% agreement). It works, but every new
evidence type (HMMs, mosaics, Helitron ends, host genes) needs new hand-placed rules, and
conflicts are settled by rule order rather than by weight of evidence.

The plan is to keep every evidence module and replace the cascade with one scored decision:

1. **Feature table.** `classify` writes one row of typed features per family
   (best bits/coverage per class from each homology source, domain-implied class, structural
   hallmarks, TSD length and support, Helitron-end support, copy number, cross-species group label).
2. **Model.** A multinomial logistic regression over top-level classes plus `non_TE`, then
   superfamily within the chosen class. Each evidence type contributes a learnable weight, so the
   model stays interpretable ("TSD 4 bp TTAA: +2.1 for DNA/PiggyBac"). Fitted offline; the weights
   ship in the refdb (`build/classifier.json`), so no new Rust dependencies.
3. **Training data.** Curated Repbase/Dfam consensi (clean labels, held out by family cluster) and the
   600k RepeatClassifier-labelled families from the 559 libraries (noisy; known weak spots such as
   SINE and RC/Helitron down-weighted). Splits are by species and clade, never by family, to avoid leakage.
4. **Decision.** Posterior per class. Below a threshold the family is `unresolved`; strong TE
   hallmarks with no winning class make a `novel_candidate`; confidence is the calibrated probability.
5. **Rollout.** Run beside the rule cascade on held-out species (`--model`), compare per class,
   and switch the default only where it is better.

- [ ] M8a feature table (`features.tsv`) from `classify`
- [ ] M8b offline fit and held-out evaluation against the rule cascade
- [ ] M8c weights in the refdb; `rep classify --model`
- [ ] M8d default switch per class once it beats the cascade

## Responsibility boundary

`rep` does discovery orchestration, copy recovery and **classification/triage**.
`reputils` does single-family biological curation (boundaries, consensus refinement) and
explicitly states that it is not a classifier. The classifier therefore lives in `rep` as
`src/classify/`, mirroring `src/copy_recovery/`. It consumes `reputils` results when they exist.

## Reference data: `rep refdb`

Since 2026-06-25 Repbase is CC0 and mirrored by Dfam with a machine-readable release
manifest (334 releases, 1998 to now, with SHA-256 checksums). `rep refdb` keeps a layered, versioned database:

| Layer | Source | Current | Contents |
| --- | --- | --- | --- |
| Repbase | `dfam.org/repbase_archive` EMBL | 31.09 (2026-09-28) | 137k consensi (KW → RepeatMasker labels), 75k CDS proteins |
| Dfam curated | `Dfam-RepeatMasker.lib.gz` | 4.0 | 26k curated consensi, including rRNA/tRNA/snRNA references |
| Dfam peptides | `RepeatPeps.lib.gz` | 4.0 | 18k labelled TE proteins |
| Local | `local/records.jsonl` | ours | custom TEs (122 proteins from the old 2024 file) and our novel finds, with provenance |

`build/` merges them into `te_nucl.fa` (146k) and `te_prot.fa` (78k) after deduplicating identical
sequences and renaming ID clashes, with every record in `records.tsv`. `rep refdb fetch` is
idempotent and only downloads new releases. Novel families flow back in through
`rep refdb add --from-classify`, so recurrences in later genomes are recognised.

The old file `repbase/RepBase29_02.2024_04_30(1).fa.gz` is superseded. Its 122 custom proteins
were the only records not in public sources, and they are now in the local layer.

## Milestones

### M0 – Evidence correctness in copy recovery (done in this pass)

Copy number, divergence and TSDs from recovered copies feed classification, so these must
be right first.

- [x] Fragment merging fused neighbouring independent copies (query ranges overlapped).
- [x] Excluded loci (`ambiguous_strand`) could still be selected by ranking.
- [x] Duplicate collapse discarded rejected loci; IDs were renumbered, so references went stale.
- [x] BLASTN used megablast with E-value 10; it now uses `dc-megablast` and `1e-10` (configurable).
- [x] The 12-column BLAST path read the wrong columns; only the 13-column layout is accepted now.
- [x] The left flank was always reported as 0 and requested intervals were wrong in `copies.tsv`.
- [x] Added per-family `loci.bed` and top-level `loci.tsv`, and fingerprints now stream the genome instead of reading it into memory.
- [x] `recover-copies` CLI uses typed parsing, validates inputs, and exposes `--blast-task`, `--evalue` and `--diversity-filter`.
- [ ] Handoff work packages 1 and 3 to 7 remain (per-stage resume, partial failures, reputils adapter, docs).

### M1 – Unknown triage core: `rep classify` (done)

- [x] Parse RepeatModeler library headers (`id#Class/Family`) and select Unknowns, or all families with `--all`.
- [x] Ordered RepBase-name → `Class/Superfamily` taxonomy with unresolved accounting.
- [x] Structural detectors in pure Rust: LTR-like terminal direct repeats, TIRs, poly-A tails,
      tandem periodicity, low complexity, six-frame ORFs, terminal motifs.
- [x] Protein homology through BLASTX (DIAMOND optional) against the RepBase proteins, aggregated per class
      with query coverage and chimera detection.
- [x] Rule-based dispositions with confidence and explicit evidence strings.
- [x] Outputs: `classification.tsv`, `summary.json`, relabelled library, `novel_candidates.fa`, ORF proteins.

### M1b – Reference database and nucleotide evidence (done in this pass)

- [x] `rep refdb fetch/register/add/build/status`, with checksum verification and offline registration.
- [x] Repbase EMBL parser (it tolerates release defects) and KW → RepeatMasker label mapping.
- [x] Local layer with origin-prefixed IDs and candidate, accepted or rejected status.
- [x] `rep classify --refdb`: BLASTX on proteins plus BLASTN on consensi, with cached indexes.
- [x] `host_or_artefact` disposition for RNA genes and artefacts, protein/nucleotide cross-checks,
      and recurrence of our own `*/Unknown` families.
- [ ] Scheduled updates, e.g. a monthly `rep refdb fetch && rep refdb build` with a changelog of new records.
- [ ] Optional Dfam uncurated (DR) families for chosen clades through FamDB partitions.
- [ ] Curation workflow for the local layer: promote candidates to accepted, rename
      `*/Unknown` families to a new class name, and export a submission to Dfam.

### M2 – Copy-informed evidence (core done)

- [x] `rep classify --copies <recover-copies dir>` (automatic in the full pipeline) reads per-family loci and copies.
- [x] Copy number, full-length copies (≥80% consensus coverage) and mean identity to the consensus (age).
- [x] TSD detection at full-length copy boundaries with ±3 bp slack. The modal exact length is checked
      against a null at positions shifted 40 bp into the flanks, and the length, motif and support are reported.
- [x] Flank conservation between copies flags truncated consensi and segmental duplications.
- [x] Decisions: a TSD consistent or inconsistent with the expected superfamily range confirms or downgrades
      known classes, suggests superfamilies for novel candidates, and dispersed multicopy repeats
      with nothing else are promoted to low-confidence novel candidates.
- [x] First real genome end to end (*Metschnikowia zobellii*): fixes for tRNAs, fragment matches and conserved flanks.
- [ ] Validate on *Cloeon dipterum* (94% Unknown), *Amanita rubescens* and the moss (jobs running).
- [ ] Copy recovery speed: about 23 s per family, dominated by MAFFT; align families in parallel.
- [ ] Boundary sharpness from the coverage profile of consensus positions across copies.
- [ ] Feed `reputils curate` outputs (refined consensus, architecture hints) back in.
- [x] Copy evidence from existing RepeatMasker annotations (`classify --rm-out <.out> --genome <fa.gz>`):
      insertions from RepeatMasker IDs, consensus coverage, divergence and genome bp, with TSDs and flanks cut
      from the genome. *Cloeon dipterum*: 2.3 min against about 6 h for copy recovery. 539 of the 559
      ensembl libraries have a `.out`.

### M3 – Protein-domain architecture (core done)

- [x] Around 60 TE-role Pfam models are extracted from a full `Pfam-A.hmm`, pressed and cached;
      `hmmscan --cut_ga` runs on the consensus ORFs (`rep classify --pfam`, pipeline `--pfam`).
- [x] Domains are ordered along the coding strand. Rules cover Gypsy and Copia (integrase position),
      Pao, DIRS, Penelope, LINE (APE+RT), Maverick, Helitron, Crypton and DNA superfamilies.
- [x] Unprecedented combinations (RT+DDE, RT+PolB, RT+helicase, YR+DDE, YR+PolB, mixed
      transposase superfamilies) give `novel_domain_architecture`. Without homology this becomes a
      `novel_candidate`; with homology it is a warning.
- [x] Benchmark on 712 intact Repbase consensi from 19 classes, using domains only and no homology:
      top-level precision 100% (0 wrong class), recall 60%, 3 superfamily errors and no false novel flags.
      Pao 37/40, Gypsy 39/40, Copia 32/40, MULE 32/40, PIF-Harbinger 34/40, LINEs ~80% to class level.
- [x] Superfamily profile HMMs (`rep refdb hmms`): labelled refdb proteins clustered with MMseqs2
      (30% identity, 50% coverage), split by label, ≥5 members aligned with MAFFT, `hmmbuild`; models
      named `hmm_<rep>#<label>`. `rep classify` scans them with `hmmsearch` (`--te-hmms`, automatic with
      `--refdb` once built) on ORFs of families without a strong (≥200 bit) BLASTX hit.
- [x] Held-out benchmark (same 712 elements, domains only; training proteins ≥80% / ≥50% identical
      to any benchmark ORF removed; 1,685 / 1,559 models):

      | set              | called | class precision | exact superfamily |
      |------------------|--------|-----------------|-------------------|
      | Pfam only        | 60.1%  | 1.000           | 59.7%             |
      | + HMMs, 80% hold | 92.4%  | 0.992           | 90.7%             |
      | + HMMs, 50% hold | 88.6%  | 0.987           | 85.7%             |

      Pfam's weak classes at the 80% hold-out: Penelope 2→40/40, DIRS 8→38, Maverick 11→34, P 11→35,
      Crypton 0→18/30 (0 at the 50% hold-out: Crypton needs close relatives). Most wrong-class calls are
      biologically ambiguous (Teratorn piggyBac–herpesvirus fusions, MuDR-7_HM with an RT domain).
- [x] Production set in `/software/team301/te_refdb/hmms`: 1,712 models, 124 labels, from 77k proteins
      (30 min on 16 cores; `.hmm` only, as hmmsearch needs no pressed files).
- [x] Real-library trial (*Agrochola circellaris*, 3,751 families, `--all --host-domains` + `.out` copies:
      6 min, 1 GB): profile hits on 195 families, all already labelled. They raise confidence
      (15 calls, e.g. PIF-Harbinger, Academ, Transib low → medium) but resolve no extra Unknowns, because
      this moth's Unknowns are almost all non-coding (14 have an ORF and no protein hit). Expect gains
      in poorly sampled taxa with diverged coding elements rather than in Lepidoptera.
- [x] Host-gene flags (`--host-domains`): ORFs of families without TE protein or domain evidence are
      scanned against all of Pfam (`hmmsearch --cut_ga`); non-TE-associated domains make the family
      `host_or_artefact` (`Other/host_gene`), or a `possible_gene_capture` warning when it has TE structure.
      Water vole: 13 of 95 Unknowns are ribosomal-protein, HMG, actin, cyclophilin, MHC and olfactory
      receptor genes (processed pseudogenes); 70 s for the whole classify run.

### M4 – Novelty discovery and negative filters (started)

`rep discover --results <dir of classify --all outputs> --out <dir> [--taxa taxa.tsv]`:
- [x] Pool all families with metadata (species, taxon, labels, dispositions, hallmarks, genome bp).
- [x] Family clusters: MMseqs2 `easy-linclust` at 80% identity and 80% coverage (both strands checked).
- [x] Groups: family clusters linked through ORF protein clusters (30% identity, 50% coverage).
      Protein clusters whose labelled members span classes are not used (bridges such as captured genes).
- [x] Label propagation to Unknowns from labelled cluster or group members (purity ≥ 0.75), with
      leave-one-out accuracy. 20-library test: 95.2% (family clusters, n=3,005) and 98.0% (groups, n=7,845).
- [x] Novel groups (no labelled member) ranked by species breadth, abundance and TE hallmarks.
- [x] Full run on all 559 v3 libraries (1.70 M families, 50 min, 6 GB): leave-one-out accuracy 96.6%
      (family clusters, n=109k) / 98.5% (groups, n=256k); 41.9k Unknowns propagated; 45k multi-species
      novel groups. Top-ranked groups are TIR-bearing non-coding families shared by 20–40 noctuid moths.
      Groups that are mostly simple/satellite/host members are no longer called novel (`non_te` column);
      cross-kingdom "novel" coding groups were host genes (7tm, actin, ribosomal S2, EGF), so
      v4 classify runs with `--host-domains`.
- [ ] Negative filters: rRNA, tRNA and snRNA (Infernal and Rfam), organellar inserts, multi-copy host genes.
- [ ] Per-candidate evidence report (alignment of members, domain maps, copy landscape).

### M5 – Feedback into masking (partly done)

- [x] The full pipeline runs `classify` between copy recovery and RepeatMasker (`--refdb`, `--pfam`),
      and RepeatMasker masks with the relabelled library. `--ltr_struct` passes `-LTRStruct` to RepeatModeler.
- [ ] Report the fraction of Unknown base pairs resolved (compare RepeatMasker tables before and after).

### M6 – `rep annotate`: one-stop TE identification for any sequence (CENSOR replacement)

GIRI's CENSOR screens arbitrary sequences against Repbase and reports classified fragments.
`rep annotate` would do the same against the merged database:

- Input any FASTA: contigs, BACs, transcripts, or a single element.
- Nucleotide and translated search against `refdb`, then fragment chaining (reusing copy
  recovery's merging) and overlap resolution by score.
- Output GFF3, BED and a CENSOR-like table (fragment, class, identity, positions in the
  reference family), plus masked FASTA and per-class base-pair summaries.
- Run the structural detectors on unmatched long regions to flag possible novel elements.
- Whole-genome masking stays with RepeatMasker in the full pipeline, using the same library.

### M7 – Validation and calibration (in progress)

Test data: `/lustre/scratch122/tol/teams/blaxter/users/mb39/repeatmodeler_ensembl` holds 559 RepeatModeler
libraries (1.70 M families, **65.5% Unknown** after RepeatClassifier); 473 have Darwin genomes for copy evidence.

Findings from the first batch (v1, 48 libraries, 2026-10-12):
- Most RepeatClassifier Unknowns have no homology to Repbase, Dfam or RepeatPeps at all, so
  homology alone cannot resolve them. Copy evidence (M2), structure and domains decide the rest.
- dc-megablast at E ≤ 1e-10 missed most diverged family relationships that RepeatClassifier
  finds. `blastn -task blastn` at E ≤ 1e-5 recovers most of them, but weak hits (median ~90 bits,
  ~26% coverage) agree with RepeatClassifier's class only ~59% of the time, so they are low-confidence calls.
- Structural detectors against composition-preserving shuffles: TIR 7.5% vs 0.6%, poly-A 2.8% vs 0.3%.
- Real-data defects handled: repeated library IDs (7 libraries) and Repbase Helitrons filed under "DNA transposon".

v3 batch (all 559 libraries, 1.70 M families, `--rm-out` copy evidence, rep 6d89abb; 2026-10-13):
- 1,113,482 RepeatClassifier Unknowns: interspersed_repeat 38.6%, unresolved 38.0%, novel_candidate
  12.3%, known_te 10.8%, satellite/simple/chimeric/host 0.5%.
- Calibrated agreement with RepeatClassifier on labelled families: high 98.5% (n=77k), medium 97.1%
  (n=223k), low 84.1% (n=181k). Overall 77.4%, the shortfall mostly families we leave Unknown.
- Main disagreement: RC/Helitron → DNA/* (9.8k), all low-confidence nucleotide-only calls (e.g. KolobokE).
  Which side is right is open; a candidate rule is to keep a specific original label when our
  only evidence is weak nucleotide homology.

RC/Helitron → DNA/* investigation (2026-10-13):
- Structure cannot settle most cases: RepeatModeler Helitron consensi rarely reach the element ends
  (TC…CTRR on 3.6% of agreed Helitrons). A copy-based end test (host A|TC…CTRR|T, null 30 bp out)
  works where ends exist (Cyaniris family-577: 0.63 vs 0.00), but some real Helitrons insert into
  A-runs and break the textbook host-T rule.
- Cause: 9 Repbase/Dfam DNA consensi drove 7.4k of the 10.6k calls, and each carries an embedded
  Helitron fragment (73–91% identity to Repbase Helitrons); 73% of the disputed hits fall inside it.
- Fix: `rep refdb mosaics` (2,824 foreign segments in 2,075 of 146k references; SINE segments
  need ≥3 supporting references; shared SINE/LINE 3' ends exempt). On the 10 most affected
  libraries: RC/Helitron → DNA 899 → 443, agreement 62.6% → 65.8%, 742 labels changed (431
  Kolobok → Helitron); ~15 regressions (Helitron references carrying DNA/LTR segments).
  The remaining calls hit references the vote missed or outvoted (P-2_CydSpl 83, Fanzor1-1_HemFu 54,
  TE0839_SO2_FAM1819 26); 92% are low-confidence nucleotide-only, which M8 should weigh properly.
- Open: 3,870 RepeatModeler "SINEs" we call LINE. RepeatModeler SINE calls are often wrong, so
  these may be our better calls; check next.

Next:
- [x] v2/v3 batch → confidence bands calibrated from agreement per evidence bin
      (`target/dev/calibrate.py` prototype; move into `rep` as a calibration table).
- [ ] Hold-out benchmark: strip labels from RepeatClassifier-classified families and score recovery per class.
- [ ] Consider `nhmmer` against Dfam curated HMMs for the most diverged relationships.
- [ ] Synthetic genomes with planted elements of known structure (end-to-end fixture with real BLAST).
