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

## Classification source candidate

The unreleased 0.2 source candidate adds a complete nucleotide single-end
database and classification workflow. It does not change the released 0.1
amplicon operations.

```text
rsomics-metagenomics database build --taxonomy taxdump --library reference.fna --mapping taxids.tsv --output reference.rsdb --threads 4
rsomics-metagenomics database inspect --database reference.rsdb
rsomics-metagenomics classify --database reference.rsdb --input reads.fastq.gz --output calls.tsv --threads 4
rsomics-metagenomics report --database reference.rsdb --input calls.tsv --output report.tsv
```

Build accepts ordered local FASTA libraries, NCBI `nodes.dmp` and `names.dmp`,
and optional merged/deleted taxonomy files. Taxa come from `kraken:taxid|ID`
tags in reference IDs or an ID/taxon TSV. The native file records input hashes,
parameters and section checksums. Existing named classifier outputs are not
replaced. Inspect reports compact-cell counts, not read counts.

Classification reads single-end FASTA/FASTQ with content-detected compression
and emits five-column calls. Report consumes those calls and the database's
taxonomy and emits the standard six columns. Equal-count sibling order is
stable by taxon ID, with observed zero-count clades before absent ones. The
candidate currently rejects zero-denominator reports; the empty-report policy
is not yet qualified for release.

The selected Kraken 2.17.1 profile retains compact fingerprints, collision LCA,
spaced minimizers and the minimum two hit groups. It uses fixed build block and
subblock defaults independent of worker count; source comparisons use the same
explicit settings. Native `.rsdb` is not Kraken `.k2d`. Paired/protein modes,
reference downloads, masking, confidence/quick modes and minimizer-data reports
are outside this slice.

See the [representative performance record](benchmarks/2026-11-05-classification.md)
for the full matched build/classify/report gate and the four-worker wall-time
loss. Source-candidate evidence is not a published 0.2 release.

This crate is licensed under either of Apache License 2.0 or MIT at your
option.
