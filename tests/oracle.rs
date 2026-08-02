use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn binary() -> PathBuf {
    env!("CARGO_BIN_EXE_rsomics-metagenomics").into()
}

fn run(mut command: Command) -> Output {
    let output = command.output().expect("required executable must run");
    assert!(
        output.status.success(),
        "command failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn require_oracle() {
    let output = run({
        let mut command = Command::new("vsearch");
        command.arg("--version");
        command
    });
    let version = format!(
        "{}{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(version.contains("v2.31.0"), "unexpected oracle: {version}");
}

fn assert_equal(left: &Path, right: &Path) {
    let left_bytes = fs::read(left).unwrap();
    let right_bytes = fs::read(right).unwrap();
    assert_eq!(
        left_bytes,
        right_bytes,
        "output mismatch\nours:\n{}\noracle:\n{}",
        String::from_utf8_lossy(&left_bytes),
        String::from_utf8_lossy(&right_bytes)
    );
}

#[test]
fn randomized_full_and_prefix_dereplication_match_vsearch() {
    require_oracle();
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("families.fasta");
    fs::write(&input, prefix_families(140)).unwrap();

    for (ours_mode, oracle_mode) in [
        ("full-length", "--derep_fulllength"),
        ("prefix", "--derep_prefix"),
    ] {
        let ours = directory.path().join(format!("ours-{ours_mode}.fasta"));
        let oracle = directory.path().join(format!("oracle-{ours_mode}.fasta"));
        run({
            let mut command = Command::new(binary());
            command.args([
                "dereplicate",
                "--input",
                input.to_str().unwrap(),
                "--output",
                ours.to_str().unwrap(),
                "--profile",
                "vsearch",
                "--mode",
                ours_mode,
                "--input-abundance",
            ]);
            command
        });
        run({
            let mut command = Command::new("vsearch");
            command.args([
                oracle_mode,
                input.to_str().unwrap(),
                "--output",
                oracle.to_str().unwrap(),
                "--sizein",
                "--sizeout",
                "--quiet",
            ]);
            command
        });
        assert_equal(&ours, &oracle);
    }
}

#[test]
fn abundance_sort_filters_and_top_match_vsearch() {
    require_oracle();
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("records.fasta");
    let ours = directory.path().join("ours.fasta");
    let oracle = directory.path().join("oracle.fasta");
    fs::write(&input, independent_records(1_200)).unwrap();

    run({
        let mut command = Command::new(binary());
        command.args([
            "sort-abundance",
            "--input",
            input.to_str().unwrap(),
            "--output",
            ours.to_str().unwrap(),
            "--profile",
            "vsearch",
            "--minimum-abundance",
            "2",
            "--maximum-abundance",
            "8",
            "--top",
            "350",
        ]);
        command
    });
    run({
        let mut command = Command::new("vsearch");
        command.args([
            "--sortbysize",
            input.to_str().unwrap(),
            "--output",
            oracle.to_str().unwrap(),
            "--sizeout",
            "--minsize",
            "2",
            "--maxsize",
            "8",
            "--topn",
            "350",
            "--quiet",
        ]);
        command
    });
    assert_equal(&ours, &oracle);
}

#[test]
fn rereplication_header_edges_match_vsearch() {
    require_oracle();
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("records.fasta");
    let ours = directory.path().join("ours.fasta");
    let oracle = directory.path().join("oracle.fasta");
    fs::write(
        &input,
        b">plain description\nACGT\n>valid;size=3\nAAAA\n>bad;size=x\nCCCC\n>duplicate;size=2;size=5\nGGGG\n",
    )
    .unwrap();

    run({
        let mut command = Command::new(binary());
        command.args([
            "rereplicate",
            "--input",
            input.to_str().unwrap(),
            "--output",
            ours.to_str().unwrap(),
            "--profile",
            "vsearch",
            "--emit-size-one",
        ]);
        command
    });
    run({
        let mut command = Command::new("vsearch");
        command.args([
            "--rereplicate",
            input.to_str().unwrap(),
            "--output",
            oracle.to_str().unwrap(),
            "--sizeout",
            "--quiet",
        ]);
        command
    });
    assert_equal(&ours, &oracle);
}

#[test]
fn empty_compatibility_input_matches_vsearch() {
    require_oracle();
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("empty.fasta");
    fs::write(&input, []).unwrap();
    for operation in ["dereplicate", "sort-abundance", "rereplicate"] {
        let ours = directory.path().join(format!("ours-{operation}.fasta"));
        let oracle = directory.path().join(format!("oracle-{operation}.fasta"));
        run({
            let mut command = Command::new(binary());
            command.args([
                operation,
                "--input",
                input.to_str().unwrap(),
                "--output",
                ours.to_str().unwrap(),
                "--profile",
                "vsearch",
            ]);
            command
        });
        let oracle_operation = match operation {
            "dereplicate" => "--derep_fulllength",
            "sort-abundance" => "--sortbysize",
            "rereplicate" => "--rereplicate",
            _ => unreachable!(),
        };
        run({
            let mut command = Command::new("vsearch");
            command.args([
                oracle_operation,
                input.to_str().unwrap(),
                "--output",
                oracle.to_str().unwrap(),
                "--sizeout",
                "--quiet",
            ]);
            command
        });
        assert_equal(&ours, &oracle);
    }
}

#[test]
fn zero_and_out_of_range_abundances_fail_in_both_implementations() {
    require_oracle();
    let directory = tempfile::tempdir().unwrap();
    for (name, header) in [
        ("zero", "x;size=0"),
        ("overflow", "x;size=9223372036854775808"),
    ] {
        let input = directory.path().join(format!("{name}.fasta"));
        fs::write(&input, format!(">{header}\nACGT\n")).unwrap();
        let ours = Command::new(binary())
            .args([
                "rereplicate",
                "--input",
                input.to_str().unwrap(),
                "--profile",
                "vsearch",
            ])
            .output()
            .unwrap();
        let oracle = Command::new("vsearch")
            .args([
                "--rereplicate",
                input.to_str().unwrap(),
                "--output",
                "/dev/null",
                "--quiet",
            ])
            .output()
            .unwrap();
        assert!(!ours.status.success());
        assert!(!oracle.status.success());
    }
}

fn prefix_families(groups: usize) -> Vec<u8> {
    let mut output = Vec::new();
    let mut state = 0x5eed_f00d_d15c_a11eu64;
    for group in 0..groups {
        let length = 72 + group % 37;
        let sequence = random_sequence(&mut state, length);
        let prefix = &sequence[..32 + group % 31];
        let lower = sequence
            .iter()
            .map(|byte| {
                if *byte == b'T' && group % 3 == 0 {
                    b'u'
                } else {
                    byte.to_ascii_lowercase()
                }
            })
            .collect::<Vec<_>>();
        extend_record(
            &mut output,
            format!("g{group:04}-prefix;size={}", group % 7 + 1).as_bytes(),
            prefix,
        );
        extend_record(
            &mut output,
            format!("g{group:04}-full;size={}", group % 11 + 1).as_bytes(),
            &sequence,
        );
        extend_record(
            &mut output,
            format!("g{group:04}-duplicate;size={}", group % 5 + 1).as_bytes(),
            &lower,
        );
    }
    output
}

fn independent_records(count: usize) -> Vec<u8> {
    let mut output = Vec::new();
    let mut state = 0x7a11_5eed_badc_0ffe;
    for index in 0..count {
        let sequence = random_sequence(&mut state, 32 + index % 83);
        extend_record(
            &mut output,
            format!("record-{index:05};size={}", index % 13 + 1).as_bytes(),
            &sequence,
        );
    }
    output
}

fn random_sequence(state: &mut u64, length: usize) -> Vec<u8> {
    (0..length)
        .map(|_| {
            *state =
                state.wrapping_add(0x9e37_79b9_7f4a_7c15).rotate_left(17) ^ 0xbf58_476d_1ce4_e5b9;
            b"ACGT"[(*state as usize) & 3]
        })
        .collect()
}

fn extend_record(output: &mut Vec<u8>, id: &[u8], sequence: &[u8]) {
    output.push(b'>');
    output.extend_from_slice(id);
    output.push(b'\n');
    output.extend_from_slice(sequence);
    output.push(b'\n');
}
