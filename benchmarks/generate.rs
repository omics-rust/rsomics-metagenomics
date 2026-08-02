use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

fn main() {
    let output = std::env::args().nth(1).expect("output directory");
    let output = Path::new(&output);
    fs::create_dir_all(output).unwrap();
    write_derep(&output.join("derep-high.fasta"), 300_000, 10_000);
    write_derep(&output.join("derep-low.fasta"), 300_000, 300_000);
    write_prefix(&output.join("prefix.fasta"), 100_000);
    write_sort(&output.join("sort.fasta"), 500_000);
    write_rereplicate(&output.join("rereplicate.fasta"), 50_000);
}

fn write_derep(path: &Path, records: usize, unique: usize) {
    let mut output = BufWriter::new(File::create(path).unwrap());
    for index in 0..records {
        let sequence = sequence(index % unique, 150);
        writeln!(output, ">record-{index:08};size={}", index % 7 + 1).unwrap();
        output.write_all(&sequence).unwrap();
        output.write_all(b"\n").unwrap();
    }
}

fn write_prefix(path: &Path, records: usize) {
    let mut output = BufWriter::new(File::create(path).unwrap());
    for index in 0..records {
        let family = index / 3;
        let length = match index % 3 {
            0 => 80,
            1 => 120,
            _ => 150,
        };
        let sequence = sequence(family, 150);
        writeln!(output, ">family-{family:08}-{length};size={}", index % 9 + 1).unwrap();
        output.write_all(&sequence[..length]).unwrap();
        output.write_all(b"\n").unwrap();
    }
}

fn write_sort(path: &Path, records: usize) {
    let mut output = BufWriter::new(File::create(path).unwrap());
    for index in 0..records {
        writeln!(output, ">sort-{index:08};size={}", index % 97 + 1).unwrap();
        output.write_all(&sequence(index, 100)).unwrap();
        output.write_all(b"\n").unwrap();
    }
}

fn write_rereplicate(path: &Path, records: usize) {
    let mut output = BufWriter::new(File::create(path).unwrap());
    for index in 0..records {
        writeln!(output, ">amplicon-{index:08};size=20").unwrap();
        output.write_all(&sequence(index, 100)).unwrap();
        output.write_all(b"\n").unwrap();
    }
}

fn sequence(index: usize, length: usize) -> Vec<u8> {
    let mut sequence = vec![b'A'; length];
    let mut identity = index;
    for byte in sequence.iter_mut().take(16) {
        *byte = b"ACGT"[identity & 3];
        identity >>= 2;
    }
    let mut state = index as u64 ^ 0x5eed_f00d_1234_5678;
    for byte in &mut sequence[16..] {
        state = state
            .wrapping_add(0x9e37_79b9_7f4a_7c15)
            .rotate_left(23)
            ^ 0xbf58_476d_1ce4_e5b9;
        *byte = b"ACGT"[(state as usize) & 3];
    }
    sequence
}
