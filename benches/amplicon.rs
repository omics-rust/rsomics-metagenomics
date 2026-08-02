use std::hint::black_box;
use std::time::Instant;

use rsomics_metagenomics::amplicon::{
    DereplicateOptions, SortOptions, dereplicate, sort_by_abundance,
};
use rsomics_seqio::OwnedRecord;

fn main() {
    let count = std::env::var("RSOMICS_BENCH_RECORDS")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(100_000usize);
    let records = fixture(count);

    let started = Instant::now();
    let dereplicated = dereplicate(
        black_box(records.clone()),
        &DereplicateOptions {
            minimum_length: 1,
            input_abundance: true,
            ..DereplicateOptions::default()
        },
    )
    .unwrap();
    let dereplicate_seconds = started.elapsed().as_secs_f64();

    let started = Instant::now();
    let sorted = sort_by_abundance(black_box(records), &SortOptions::default()).unwrap();
    let sort_seconds = started.elapsed().as_secs_f64();

    println!(
        "records={count} clusters={} dereplicate_seconds={dereplicate_seconds:.6} sort_seconds={sort_seconds:.6} sorted={}",
        dereplicated.records().len(),
        sorted.records().len()
    );
}

fn fixture(count: usize) -> Vec<OwnedRecord> {
    (0..count)
        .map(|index| {
            let unique = index % 10_000;
            let mut sequence = vec![b'A'; 250];
            for (offset, byte) in unique.to_le_bytes().into_iter().enumerate() {
                sequence[offset] = b"ACGT"[usize::from(byte & 3)];
                sequence[offset + 8] = b"ACGT"[usize::from((byte >> 2) & 3)];
            }
            OwnedRecord {
                id: format!("record-{index:08};size={}", index % 7 + 1).into_bytes(),
                seq: sequence,
                qual: None,
            }
        })
        .collect()
}
