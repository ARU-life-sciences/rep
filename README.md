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

codex resume 019fc7db-f0d4-75a0-ba31-5ae6c5133d38

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
