# rsomics-metagenomics

`rsomics-metagenomics` provides one coherent CLI and Rust library for
abundance-aware amplicon workflows. The initial release implements full-length
and prefix dereplication, abundance sorting, and checked rereplication.

```text
rsomics-metagenomics dereplicate --input reads.fasta --output uniques.fasta
rsomics-metagenomics sort-abundance --input uniques.fasta --output sorted.fasta
rsomics-metagenomics rereplicate --input sorted.fasta --output reads.fasta
```

Strict mode is the default. It requires one valid, non-zero `size` attribute
where abundance is consumed, rejects duplicate and malformed values, checks
all arithmetic, and preserves full identifiers. `--profile vsearch` provides
the selected VSEARCH 2.31.0 label, normalization, missing-abundance, and output
behavior used by compatibility tests.

Named outputs are transactional. Rereplication computes the complete output
record and byte counts before it opens the transactional output and refuses
the operation when either configured expansion limit would be exceeded.

This crate is licensed under either of Apache License 2.0 or MIT at your
option.
