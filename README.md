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

The Rust API also exposes `CurationConfig` and `run_curation_with_config` to customize flank length, maximum hit distance, orientation threshold, and hits per query.
