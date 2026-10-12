# rep roadmap: exhaustive resolution of Unknown TE families and discovery of novel classes

## Goal

1. **Every** RepeatModeler/RepeatMasker family labelled `Unknown` receives an
   evidence-backed *disposition*. No family stays unexplained without a recorded reason.
2. Families that look like real transposable elements but match no known class are
   surfaced and ranked as **candidate novel TE classes**.

A disposition is one of:

| Disposition | Meaning |
| --- | --- |
| `known_te` | Protein (and later nucleotide) homology to a classified TE; relabelled `Class/Superfamily` |
| `chimeric` | Homology to two or more incompatible TE classes in different regions |
| `simple_repeat` / `satellite` | Tandem periodicity explains the family (period recorded) |
| `low_complexity` | Compositionally biased, non-periodic |
| `novel_candidate` | TE hallmarks (terminal repeats, coding capacity, poly-A, TSDs, copy number) without known homology |
| `host_or_artefact` | Multi-copy host gene, segmental duplication, organellar or RNA-gene derived (milestone 4) |
| `unresolved` | No evidence either way; the missing evidence is recorded |

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

### M1 – Unknown triage core: `rep classify` (started)

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

### M2 – Copy-informed evidence

- Run copy recovery on Unknowns and read `copies.tsv` and the aligned copies.
- Copy number and mean divergence from the consensus (age).
- TSD detection from the flanks of selected copies (length mode and motif): 2 bp for Tc1/mariner (TA),
  4 to 6 bp for LTR/Gypsy and Copia, 8 bp for hAT, 9 to 11 bp for MULE, none for Helitrons, and variable for LINEs.
- Boundary sharpness, meaning where the copy alignment loses coverage, which separates TEs from segmental duplications.
- Feed `reputils curate` outputs (refined consensus, architecture hints) back in where they exist.

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
- [ ] Poorly detected by Pfam: Crypton (0/30), Penelope (2/40), DIRS (8/40), Maverick (11/40), P (11/40).
      Build our own superfamily-specific profile HMMs from the 78k labelled refdb proteins
      (cluster → align → `hmmbuild`) and scan them alongside the Pfam subset.
- [ ] Host-gene flags: scan ORFs that lack TE domains against full Pfam (kinases, zinc fingers, ...).

### M4 – Novelty discovery and negative filters

- All-vs-all clustering of residual Unknowns (BLASTN or MMseqs2) into Unknown groups, so
  the same novel element found many times by RepeatModeler is treated as one class candidate.
- Protein-level clustering of ORFs across groups, to find novel element families that share a protein.
- Novel domain combinations, for example YR together with a helicase, or a transposase next to an unusual accessory domain.
- Novelty score: hallmarks × copy number × (1 − best homology) × structural consistency.
- Negative filters: rRNA, tRNA and snRNA (Infernal and Rfam), organellar inserts, segmental
  duplications, and multi-copy host genes.
- Per-candidate HTML or TSV evidence report.

### M5 – Feedback into masking

- Relabelled library into RepeatMasker for a second round, reporting the fraction of
  Unknown base pairs resolved.
- Route the full pipeline through `classify` between copy recovery and RepeatMasker.

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

### M7 – Validation

- Hold out families RepeatModeler *did* classify, strip their labels, and measure how well
  dispositions recover them (precision and recall per class).
- Use the 122 custom labelled proteins as positive controls for the taxonomy and homology path.
- Synthetic genomes with planted elements of known structure, including an end-to-end fixture with real BLAST.
