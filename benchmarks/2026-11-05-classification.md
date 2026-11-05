# Native compact classification performance

Apple M2, eight cores, 8 GiB RAM, Darwin27.2 arm64, Rust1.91. Exact upstream
Kraken2.17.1 source `5e2aa928d00b96d61f204d517437637863da1d8c` is built
with Apple Clang and OpenMP. No builds or other task benchmarks ran during
timing. This is one native host, not a cross-platform speed claim.

The controlled workload uses 55,946,218 actual reference bases: the full
1000 Genomes human b37 chromosome22 followed by NCBI ASM584v2 E. coli K-12
MG1655. Its42-node taxonomy subset is derived from actual NCBI lineage records.
One million 150-base queries use fixed seed `0x52534F4D49435301`:450,000 valid
human windows,450,000 E. coli windows,80,000 random controls and20,000 all-N.
All900,000 source-derived coordinates were verified before timing. This is a
synthetic query workload on real references, not biological accuracy truth or
a benchmark of loading the complete NCBI taxonomy.

Both builds use k35/l31/spaces7, sampled capacity, load factor0.7 and explicit
block/subblock16384. No masking or hash downsampling. The source build includes
the real wrapper's estimator and construction; both inputs have exactly the
same order. Fixed partitioning is necessary for the matched compatibility
profile, so the four-worker source lane is not Kraken's automatic subblock
default. The query lane includes native classification and the standalone
report; upstream classification uses its integrated `--report`.

One warmup and five alternating measured pairs per worker count/workflow:
24 total pairs,60 subprocesses. Wall is monotonic invocation time; native
sequential classify/report wall and CPU times are summed. RSS is the maximum
observed process RSS from Darwin `time -l`, max across native sequential
stages, not aggregate concurrent-process memory. Checks are outside timing.

| Workflow | Workers | Native wall median s | Kraken wall median s | Native / Kraken median RSS bytes | Strict advantage |
|---|---:|---:|---:|---:|---|
| build |1|4.373175|7.868703|131727360 / 175063040|wall and RSS|
| build |4|3.757125|8.039115|135905280 / 175128576|wall only|
| classify/report |1|8.473688|9.122271|75382784 / 140378112|wall and RSS|
| classify/report |4|5.387512|4.896346|75546624 / 309411840|RSS only|

Strict means every measured native observation is below every source
observation on that axis. All four lanes pass at least one axis. The four-
worker build RSS ranges overlap, including a lower source observation. Four-
worker classify/report has about10% slower median wall, with overlapping
ranges, but approximately76% lower median RSS. Native classify/report scales
about1.57 times from one to four workers with disjoint wall ranges; native
build's worker ranges overlap. Do not claim general throughput or scaling.

All full outputs match across both tools, workers and trials. The decoded
taxonomies and all16,653,531 table cells match;11,654,303 are occupied with
six value bits. Both report totals account for all one million reads. Complete
classification SHA-256:
`b4a0f985d8b645ef82bf3f89c6ab14264b38e65efb72eee0304e4d7b9be8ee5c`.
Complete report SHA-256:
`127755e83cb9f9d316cc307b2ea7dbcb30bd6deb65909f449b02dafc64655bd5`.

Measured binary SHA-256:
`42823c26093f12cc74e04e3edfaab4a532ad890aedb3060071f7771337de1f36`.
It is the working implementation subsequently committed at `9bfb179`, before
the 0.2 version change. Every production source body was independently matched
to the measurement receipt. No claim is made that a later package binary was
timed.

Original packet:
`/Volumes/KIOXIA/Developments/tmp/metagenomics-classify-perf-20261105.HaqAWrwA`.
The `preparation` manifest and verification freeze the inputs. `campaign`
contains all argv, outputs, individual wall/CPU/RSS measurements, original
database files and the retained measured executable.

- input `MANIFEST.json`: `52c34860fa0a1f7d42ec2fa23a204ca99bcaa1aece1ac66ce4434520f67d190f`;
- `IDENTITY.json`: `baa808418f5d045dabb7042c7caf4d522b4bf15f9f4ab69d53f0d3c6658d40d0`;
- `samples.json`: `776857c6d8d28bc2cb83e9e073f72b10ea19e8ba80715d569a016e879acfa372`;
- `COMPLETE.json`: `570bd9b147e84bb351fe5c7ba838626ee7a2dab9ab20dab88fa84504cd40d0ee`.
