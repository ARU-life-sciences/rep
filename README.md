# Automate a repeat annotation pipeline

This tool runs RepeatModeler, curates each discovered repeat family against the genome, and then runs RepeatMasker.

The curation stage:

1. Builds a nucleotide BLAST database from the staged genome.
2. Searches the RepeatModeler `consensi.fa.classified` library against it.
3. Keeps the 20 best hits per family, merges nearby fragments, and adds 2 kb flanks.
4. Reverse-complements regions whose hits are predominantly on the minus strand.
5. Writes per-family FASTA files and MAFFT alignments under `curation/`.

Required executables are resolved from `PATH`. Custom installations can be selected with:

- `REP_REPEATMASKER`
- `REP_REPEATMODELER`
- `REP_BUILD_DATABASE`
- `REP_MAKEBLASTDB`
- `REP_BLASTN`
- `REP_MAFFT`

Use `--curation_target <FASTA>` to search the RepeatModeler consensus against a different genome; gzip-compressed targets are decompressed into the curation workspace. Use `--curation_only` to resume curation from an existing RepeatModeler output without rerunning discovery.

The Rust API also exposes `CurationConfig` and `run_curation_with_config` to customize flank length, maximum hit distance, orientation threshold, and hits per query.

## Copy recovery

`rep recover-copies --queries consensi.fa.classified --genome genome.fa --out copies/`
performs BLASTN search, coherent fragment merging, independent-locus selection, clipped
sequence extraction, strand orientation, and MAFFT alignment. Coordinates in BLAST evidence
are 1-based inclusive; extraction metadata uses 0-based half-open intervals. Minus-strand
copies are reverse-complemented into query-relative FASTA orientation while genomic strand
and coordinates remain in metadata.

Each family directory contains `query.fa`, `hits.tsv`, `loci.tsv`, `copies.tsv`,
`copies.raw.fa`, `copies.extended.fa`, `copies.aligned.fa`, and `manifest.json`.
Families with no selected loci receive an empty alignment (`no_usable_copies`); one-copy
families are copied through as `single_copy`; multi-copy families are aligned with MAFFT.
Manifests use the versioned `rep.copy-recovery.v1` contract and include stage fields.
The legacy `CurationConfig` API remains available for compatibility. Boundary inference,
structural interpretation, and protein-domain analysis remain downstream responsibilities
of `reputils`.

## Classifying Unknown families

`rep classify --library consensi.fa.classified --proteins RepBase_proteins.fa.gz --out classify/`
gives every `Unknown` family an evidence-backed disposition: `known_te`, `chimeric`,
`simple_repeat`, `satellite`, `low_complexity`, `novel_candidate` (TE hallmarks such as
terminal direct/inverted repeats, poly-A tails or coding capacity without known homology)
or `unresolved` (with the missing evidence recorded). Evidence comes from BLASTX (or
`--backend diamond`) against TE proteins, whose names are mapped to RepeatMasker classes,
and from built-in structural detectors. Use `--all` to also classify labelled families and
report agreement with their existing labels.

Outputs: `classification.tsv`, `families.json`, `summary.json`, `library.reclassified.fa`
(Unknowns relabelled for RepeatMasker), `novel_candidates.fa`, `orfs.faa` and
`unclassified_subjects.tsv` (protein hits whose names the taxonomy could not place).
The `host_or_artefact` disposition marks families matching RNA genes or other non-TE references.
Executables can be overridden with `REP_BLASTX`, `REP_MAKEBLASTDB` and `REP_DIAMOND`.

## Reference database

`rep refdb` maintains a merged, updateable TE reference built from Repbase (CC0 since June
2026, mirrored by Dfam), Dfam curated families and peptides, and a local layer of our own families:

```bash
rep refdb fetch  --db te_refdb            # latest Repbase + Dfam; only new releases are downloaded
rep refdb add    --db te_refdb --fasta custom_te.faa --kind protein --origin lab --status accepted
rep refdb build  --db te_refdb            # build/te_nucl.fa, build/te_prot.fa, build/records.tsv
rep classify --library consensi.fa.classified --refdb te_refdb --out classify/
rep refdb add    --db te_refdb --from-classify classify/ --origin <genome> --status candidate
```

With `--refdb`, classification uses both protein (BLASTX) and nucleotide (BLASTN) homology.
Search indexes are cached in `build/index/` until the next build. `rep refdb register`
registers already-downloaded releases for offline use.

See `ROADMAP.md` for the plan towards exhaustive Unknown resolution and novel-class discovery.
