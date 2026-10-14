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

- [x] M8a feature table (`features.tsv`, schema `rep.features.v1`) from every `classify` run: 105 columns.
      Targets and current call (original label, label, disposition, confidence); per class
      (DNA, RC, LTR, LINE, SINE, PLE, Retroposon) BLASTX bits/coverage, profile-HMM bits, BLASTN
      bits/coverage/identity (plus nonTE); mosaic and chimera flags; domain-implied class and 15 domain
      roles; host-gene domains; structure (length, GC, LTR, TIR, tails, tandem, ORFs, termini); copies
      (loci, full-length, identity, genome bp, TSD length/support/null, flank similarity, Helitron ends).
- [x] Training data, libraries: v4 batch over the 559 libraries (rep c34fbe6; `--rm-out`, `--host-domains`,
      mosaics, HMMs): 1.70M feature rows, 585k with a RepeatClassifier label (LINE 241k, DNA 156k,
      LTR 85k, RC 63k, SINE 34k, PLE 3k); 121 CPU-hours.
- [x] Training data, curated: refdb consensi clustered (≥80% identity, member ≥80% covered), clusters
      in 5 folds (~29k consensi each); fold k's reference drops its consensi and the proteins they
      encode (≥80% identity over ≥50%; ~22% of proteins), with HMMs rebuilt. 139k labelled rows.
      The cascade on them: calls 80.8%, agreement 0.986; LTR recall only 0.715 (solo `-LTR` entries
      without hallmarks), SINE precision 0.971 (0.854 against RepeatModeler's SINE labels).
- [x] M8b offline fit and held-out evaluation against the rule cascade: `scripts/classifier/fit.py`
      (numpy/scipy multinomial logistic regression; log1p scores and counts, missing indicators,
      standardised; class weights ∝ frequency^-0.5; trained on library and curated rows).
      Libraries held out by order (5 folds; Lepidoptera alone is one): at the cascade's 83.8% coverage,
      agreement 0.941 vs 0.927; at the cascade's agreement, coverage 86.7%. Curated consensi held
      out by cluster fold: 0.998 vs 0.986 at equal coverage; 96.2% vs 80.8% coverage at equal agreement.
      Gate: the cascade's non-TE calls (simple repeats, satellites, RNA and host genes) override the
      model, restoring nonTE recall 0.51 → 0.97 at no cost elsewhere. The model removes most of the
      cascade's systematic confusions against RepeatModeler labels (RC→DNA 6,750 → 2,291,
      SINE→LINE 3,976 → 1,956) and adds LINE→SINE (823 → 2,046) and LINE→RC (659 → 1,573).
      Superfamily: given the right class, the best-hit superfamily is right in 96.4% of curated
      consensi (errors within clades: CR1-clade → L2, ERV1/ERVL/ERVK), so the model chooses the class
      and the superfamily comes from the best hit within it. Posteriors are calibrated as they are
      (held-out agreement 0.70–0.80: 78%, 0.90–0.95: 94%, ≥0.99: 98.2%; the shortfall at the top is label
      noise). Gated, on curated consensi: 0.997 at the cascade's coverage, 96.7% coverage at its agreement.
      Final fit on all 724k rows: `build/classifier.json` (v2, 2026-10-13).
- [x] M8c weights in the refdb (`build/classifier.json`, shown by `rep refdb status`); `classify` loads it
      automatically with `--refdb` (or `--model <json>`) and adds model_class, model_label,
      model_probability and model_confidence (high ≥0.99, medium ≥0.9, low ≥0.7, else none) to
      `classification.tsv`, with agreement against the cascade in `summary.json`. The cascade's decision
      is unchanged. Rust and Python inference agree to 7e-5 (TSV rounding).
- [x] M8d the model decides where it beats the cascade (`scripts/classifier/switch_report.py`): its call
      replaces the cascade's at ≥0.7 where the cascade made a TE call and at ≥0.9 (medium) where it
      made none, except that the cascade keeps its non-TE calls, its calls backed by a blastx hit in
      their own class, and families the model calls Retroposon or nonTE. Held out: libraries
      83.8%/0.927 → 84.0%/0.936 (calls/agreement), curated consensi 80.8%/0.986 → 95.0%/0.989; LTR
      recall on curated 0.715 → 0.957, RC 0.61 → 0.68 and SINE 0.57 → 0.63 on libraries for about one
      point of precision; Retroposon and nonTE unchanged (an unrestricted switch at 0.7 cost them 16
      and 3–7 points of precision). Over the 559 Ensembl libraries (v5) the switch mostly
      relabels weak nucleotide-only cascade calls (DNA→RC 3,691, LINE→SINE 1,433, LTR→RC 830), where a
      30-family review sides with the RepeatModeler label (which the model never sees) 21–28 times;
      the model's LINE→RC calls were wrong where the cascade had a LINE protein hit, which the rule
      now protects. It resolves few Unknowns (1.4k at ≥0.9 of 985k the cascade leaves): those need
      new evidence, not a new decision rule. `classification.tsv` gains `decided_by` (the cascade's
      label stays in the evidence, `model:<label>@p;cascade:<label>`); features.tsv still records the
      cascade's call, so training data is unaffected.

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
- [x] BLASTN used megablast with E-value 10; it now uses `blastn` at `1e-5` (configurable; `dc-megablast`
      at `1e-10` until 2026-10-13, see M2).
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
- [x] Copy-recovery search on *Cloeon dipterum* (898 families): `blastn` at 1e-5 against `dc-megablast` at 1e-10
      finds a median 115.5 vs 58.5 loci and 29 vs 21 full-length copies per family (Spearman 0.86/0.93),
      identity 81.7% vs 83.6%; same TSD length in 71 families, 17 called by one only; 7.1 h vs 5.4 h.
      Against the RepeatMasker `.out` (67 loci, 15 full-length) it evaluates TSDs, flanks and Helitron
      ends in more families. Default switched to `blastn` at 1e-5.
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
- [x] Unknowns from related families (`classify` second pass, `rep relatives`). The 984k Unknowns that
      the cascade and model leave uncalled (v5) have essentially no refdb homology (1% any nucleotide hit,
      no protein or domain hits), are short (76% < 500 bp, 0.6% with an ORF ≥ 200 aa) and real
      (81% with ≥ 5 full-length copies; 62 Gb of genome). Many are pieces of elements that RepeatModeler
      split, or relatives of families called in related genomes. `classify` now searches uncalled Unknowns
      (blastn) against the families the run called as TEs and against `--relatives <taxon>.fa`, which
      `rep relatives --results <dir> --taxa <tsv> --rank order` writes from many runs (references: cascade
      TE calls and model decisions, never relative-labelled families). The best match's label is taken,
      medium confidence at ≥ 60% coverage, low below; `decided_by = relative`. Held out on labelled
      families (40 libraries), agreement 0.93–0.99 at ≥ 60% coverage and 0.81–0.95 below (LINE/SINE and
      RC/DNA confusions in partial matches). In Rust on the same libraries with per-order relatives from
      v5 (605k references, 53 orders): 55.8% of all Unknowns resolved (29.2k medium, 15.6k low;
      DNA 20.1k, LINE 11.6k, RC 5.3k, LTR 4.1k, SINE 3.1k), at a median 669 s per library against 617 s.
- [x] Long-monomer satellites from copy arrangement: an uncalled family whose RepeatMasker insertions
      average ≥ 3× its consensus length (≥ 10 loci) is a `satellite` (`tandem_array:<ratio>`). Its consensus
      holds one or two monomers of 300–1,000 bp, so it shows no periodicity (the largest open family,
      *Bombylius discolor* rnd-6_family-1284, was "LTR/Unknown": two copies of an 803 bp unit, 186 loci over
      24.7 Mb). High-confidence known TEs: ratio 99th percentile 1.01, 0.07% above 1.5. On the 40 libraries:
      48 open Unknowns (76 Mb, 5% of open bp) and 2 of ~4,100 open RepeatModeler-labelled TE families.
- [x] v6 batch (all 559 libraries; model decisions, relatives pass, tandem arrays): 1.70 M families,
      decided by cascade 1.07 M, relatives 610k, model 15.7k; 30.3% of the 1.11 M RepeatModeler Unknowns
      remain Unknown; 1,262 tandem-array satellites (2.5 Gb), 212 of them RepeatModeler "TEs" (median
      92 bp, no homology). `discover` on v6: 5.0k multi-species groups where every member is still open
      (12.6k families, 375 Mb; mostly two-species insect pairs, none with TIRs: the TIR-rich groups are now
      labelled DNA by relatives).
- [x] Cluster consensus over relatives. Relative labels agreed with the confident labels of their
      discover group only 0.74 of the time on v6 (0.87 for matches within the library), against 0.985
      leave-one-out for group propagation. References are now our confident calls (model decisions,
      cascade calls above low), else RepeatClassifier's label where we did not contradict its class, never
      relative labels; members whose propagated label differs in class go to `relabels.tsv`, which
      `rep relabel --overrides` applies.
- [ ] Negative filters: organellar inserts, multi-copy host genes. Rfam (cmscan, all 4,227 models,
      E ≤ 1e-5) does not help: 0.6% of RepeatModeler SINE/tRNA families hit a tRNA model (heads too
      degraded), and uncalled Unknowns hit no tRNA, 5S, 7SL, rRNA or snRNA model (0.3% hit plant miRNA
      precursors, themselves often MITE-derived).
- [ ] Per-candidate evidence report (alignment of members, domain maps, copy landscape).

### M5 – Feedback into masking (partly done)

- [x] The full pipeline with `--refdb` follows the batch route: RepeatModeler → RepeatMasker with its
      library → `classify --all` (copy evidence from the `.out`, `--host-domains` with `--pfam`,
      `--relatives`) → `relabel` into `results/annotation/`. Copy-recovery curation runs with `--curate`
      (always without `--refdb`). A finished RepeatModeler library and RepeatMasker annotation are reused,
      and `--rmo_only` stops after RepeatModeler, so a second pass with `--relatives` reruns only classify.
- [x] `rep relabel --classify <dir> --rm-out <.out> [--genome-sizes] [--overrides relabels.tsv --sample]`
      rewrites the class/family column of an existing annotation (no remasking) and reports base pairs per
      class before and after (union of intervals), Unknown bp resolved by class and by deciding stage.
      *Bombylius discolor*: Unknown 23.2% → 4.3% of the genome, Satellite 0.25% → 12.8%, SINE 0.17% → 0.01%.
- [ ] Batch numbers over the 559 libraries, with and without the discover relabels.

### Mosses (started 2026-10-14)

Bryophytes have no RepeatModeler libraries in the ToL or Ensembl data (7 ToL species on the Ensembl FTP),
so the panel runs the full pipeline: 35 genomes (`rep_moss_panel/panel.tsv`), one or more per moss order
plus two liverworts, 15.4 Gb in total. Stage 1 `--rmo_only` (RepeatModeler -LTRStruct), stage 2 the
`--refdb` pipeline, then `rep relatives` over the panel and a classify second pass with `--relatives`,
then `discover`. The earlier *Ptychomitrium polyphyllum* run (rep 96af2f3) left 387 of 430 families
Unknown with only 4 known TEs.

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
