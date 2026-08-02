use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};

fn binary() -> PathBuf {
    env!("CARGO_BIN_EXE_rsomics-metagenomics").into()
}

fn run(args: &[&str]) -> Output {
    Command::new(binary()).args(args).output().unwrap()
}

#[test]
fn strict_sort_requires_abundance_and_preserves_destination() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.fasta");
    let output = directory.path().join("output.fasta");
    fs::write(&input, b">missing\nACGT\n").unwrap();
    fs::write(&output, b"sentinel\n").unwrap();
    let result = run(&[
        "sort-abundance",
        "--input",
        input.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);
    assert_eq!(result.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&result.stderr).contains("has no size abundance"));
    assert_eq!(fs::read(&output).unwrap(), b"sentinel\n");
}

#[test]
fn strict_dereplication_rejects_malformed_ignored_abundance() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.fasta");
    let output = directory.path().join("output.fasta");
    fs::write(
        &input,
        b">good\nACGTACGTACGTACGTACGTACGTACGTACGT\n>bad;size=no\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
    )
    .unwrap();
    fs::write(&output, b"sentinel\n").unwrap();
    let result = run(&[
        "dereplicate",
        "--input",
        input.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);
    assert_eq!(result.status.code(), Some(1));
    assert_eq!(fs::read(&output).unwrap(), b"sentinel\n");
}

#[test]
fn rereplication_budget_is_enforced_before_transaction() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.fasta");
    let output = directory.path().join("output.fasta");
    fs::write(&input, b">many;size=50\nACGT\n").unwrap();
    fs::write(&output, b"sentinel\n").unwrap();
    let result = run(&[
        "rereplicate",
        "--input",
        input.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
        "--maximum-records",
        "49",
    ]);
    assert_eq!(result.status.code(), Some(2));
    assert_eq!(fs::read(&output).unwrap(), b"sentinel\n");
}

#[test]
fn named_json_output_uses_shared_envelope() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.fasta");
    let output = directory.path().join("output.fasta");
    fs::write(&input, b">one;size=2\nACGT\n").unwrap();
    let result = run(&[
        "sort-abundance",
        "--input",
        input.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
        "--json",
    ]);
    assert!(result.status.success());
    let report: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(report["tool"], "rsomics-metagenomics");
    assert_eq!(report["result"]["operation"], "sort-abundance");
    assert_eq!(report["result"]["output_records"], 1);
}

#[test]
fn output_alias_is_rejected() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.fasta");
    fs::write(&input, b">one;size=2\nACGT\n").unwrap();
    let result = run(&[
        "sort-abundance",
        "--input",
        input.to_str().unwrap(),
        "--output",
        input.to_str().unwrap(),
    ]);
    assert_eq!(result.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&result.stderr).contains("also an input path"));
}
