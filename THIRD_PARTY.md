# Behavior provenance

The `vsearch` compatibility profile is tested against VSEARCH 2.31.0 from the
official Torognes VSEARCH release. VSEARCH is available under the BSD 2-Clause
or GPL-3.0 license; this project selects the BSD-compatible behavior reference.

The consolidated Rust implementation and its historical rsomics source assets
are team-owned. VSEARCH names identify compatibility targets and do not make
VSEARCH a runtime dependency.

The native classifier's minimizer, compact hash, taxonomy assignment and
standard report behavior derive from Kraken 2.17.1, exact source revision
`5e2aa928d00b96d61f204d517437637863da1d8c` of
[DerrickWood/kraken2](https://github.com/DerrickWood/kraken2). See the original
[Kraken 2 MIT license](LICENSES/kraken2-MIT). Kraken 2 is a test oracle, not a
runtime dependency. The native container does not claim `.k2d` interchange.
