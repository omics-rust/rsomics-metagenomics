use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/classification")
        .join(name)
}

fn tool() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-metagenomics"))
}

fn succeeds(command: &mut Command) -> Output {
    let result = command.output().unwrap();
    assert!(
        result.status.success(),
        "{command:?}\n{}",
        String::from_utf8_lossy(&result.stderr)
    );
    result
}

#[test]
fn database_classify_report_matches_pinned_kraken() {
    let work = tempfile::tempdir().unwrap();
    for (label, extra, input, calls, report) in [
        (
            "basic",
            vec![],
            "reads.fastq",
            "calls-basic.tsv",
            "report-basic.tsv",
        ),
        (
            "dense",
            vec![
                "--capacity",
                "1024",
                "--block-size",
                "256",
                "--subblock-size",
                "64",
                "--minimum-bits-for-taxid",
                "3",
            ],
            "unknown.fastq",
            "calls-v3.tsv",
            "report-v3.tsv",
        ),
        (
            "collision",
            vec![
                "--capacity",
                "1024",
                "--block-size",
                "256",
                "--subblock-size",
                "64",
                "--minimum-bits-for-taxid",
                "28",
            ],
            "unknown.fastq",
            "calls-v28.tsv",
            "report-v28.tsv",
        ),
    ] {
        for workers in ["1", "4"] {
            let db = work.path().join(format!("{label}-{workers}.rsdb"));
            let output = work.path().join(format!("{label}-{workers}.calls"));
            let summary = work.path().join(format!("{label}-{workers}.report"));
            succeeds(
                tool()
                    .args(["database", "build", "--taxonomy"])
                    .arg(fixture("taxonomy"))
                    .arg("--library")
                    .arg(fixture("reference.fna"))
                    .arg("--output")
                    .arg(&db)
                    .args(["--threads", workers])
                    .args(&extra),
            );
            succeeds(
                tool()
                    .args(["classify", "--database"])
                    .arg(&db)
                    .arg("--input")
                    .arg(fixture(input))
                    .arg("--output")
                    .arg(&output)
                    .args(["--threads", workers]),
            );
            assert_eq!(
                std::fs::read(&output).unwrap(),
                std::fs::read(fixture(calls)).unwrap(),
                "{label}/{workers}"
            );
            succeeds(
                tool()
                    .args(["report", "--database"])
                    .arg(&db)
                    .arg("--input")
                    .arg(&output)
                    .arg("--output")
                    .arg(&summary)
                    .args(if label == "basic" {
                        &["--report-zero-counts"][..]
                    } else {
                        &[]
                    }),
            );
            assert_eq!(
                std::fs::read(summary).unwrap(),
                std::fs::read(fixture(report)).unwrap(),
                "{label}/{workers}"
            );
            let metadata = succeeds(
                tool()
                    .args(["--json", "database", "inspect", "--database"])
                    .arg(&db),
            );
            let metadata: serde_json::Value = serde_json::from_slice(&metadata.stdout).unwrap();
            assert_eq!(
                metadata["result"]["capacity"],
                if label == "basic" { 13165 } else { 1024 }
            );
            assert_eq!(
                metadata["result"]["occupied"],
                match label {
                    "basic" => 700,
                    "dense" => 702,
                    _ => 662,
                }
            );
        }
    }
}

#[test]
fn invalid_database_and_late_query_do_not_commit_outputs() {
    let work = tempfile::tempdir().unwrap();
    let bad_taxonomy = work.path().join("invalid-taxonomy");
    std::fs::create_dir(&bad_taxonomy).unwrap();
    let nodes = std::fs::read_to_string(fixture("taxonomy/nodes.dmp"))
        .unwrap()
        .replace("10\t|\t2\t|", "10\t|\t11\t|");
    std::fs::write(bad_taxonomy.join("nodes.dmp"), nodes).unwrap();
    std::fs::copy(
        fixture("taxonomy/names.dmp"),
        bad_taxonomy.join("names.dmp"),
    )
    .unwrap();
    let rejected = work.path().join("rejected.rsdb");
    let result = tool()
        .args(["database", "build", "--taxonomy"])
        .arg(&bad_taxonomy)
        .arg("--library")
        .arg(fixture("reference.fna"))
        .arg("--output")
        .arg(&rejected)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("cycle"));
    assert!(!rejected.exists());
    let db = work.path().join("reference.rsdb");
    succeeds(
        tool()
            .args(["database", "build", "--taxonomy"])
            .arg(fixture("taxonomy"))
            .arg("--library")
            .arg(fixture("reference.fna"))
            .arg("--output")
            .arg(&db),
    );
    let input = work.path().join("malformed.fastq");
    let mut bytes = std::fs::read(fixture("reads.fastq")).unwrap();
    bytes.extend_from_slice(b"@broken\nACGT\n+\nIII\n");
    std::fs::write(&input, bytes).unwrap();
    let output = work.path().join("calls.tsv");
    let result = tool()
        .args(["classify", "--database"])
        .arg(&db)
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(!output.exists());
    let mut bytes = std::fs::read(&db).unwrap();
    *bytes.last_mut().unwrap() ^= 1;
    std::fs::write(&db, bytes).unwrap();
    std::fs::write(&output, b"previous result\n").unwrap();
    let result = tool()
        .args(["classify", "--database"])
        .arg(&db)
        .arg("--input")
        .arg(fixture("reads.fastq"))
        .arg("--output")
        .arg(&output)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("integrity"));
    assert_eq!(std::fs::read(output).unwrap(), b"previous result\n");
}

#[test]
fn live_kraken_preserves_compact_cells_and_connected_outputs() {
    let oracle = std::env::var_os("RSOMICS_KRAKEN2_DIR")
        .map(PathBuf::from)
        .or_else(|| {
            std::env::var_os("PATH").and_then(|paths| {
                std::env::split_paths(&paths).find(|path| path.join("kraken2").is_file())
            })
        })
        .expect("Kraken2 2.17.1 source oracle is required");
    let version = succeeds(Command::new(oracle.join("kraken2")).arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("Kraken version 2.17.1\n"));
    let work = tempfile::tempdir().unwrap();
    let source_db = work.path().join("source");
    std::fs::create_dir(&source_db).unwrap();
    let mapping = work.path().join("mapping.tsv");
    std::fs::write(
        &mapping,
        b"ref11|kraken:taxid|11\t11\nref12|kraken:taxid|12\t12\nref21|kraken:taxid|21\t21\n",
    )
    .unwrap();
    let library = std::fs::File::open(fixture("reference.fna")).unwrap();
    succeeds(
        Command::new(oracle.join("build_db"))
            .args([
                "-k",
                "35",
                "-l",
                "31",
                "-S",
                "1111111111111111101010101010101",
                "-c",
                "1024",
                "-p",
                "4",
                "-B",
                "256",
                "-b",
                "64",
                "-r",
                "28",
            ])
            .arg("-n")
            .arg(fixture("taxonomy"))
            .arg("-m")
            .arg(mapping)
            .arg("-H")
            .arg(source_db.join("hash.k2d"))
            .arg("-t")
            .arg(source_db.join("taxo.k2d"))
            .arg("-o")
            .arg(source_db.join("opts.k2d"))
            .stdin(library),
    );
    let native_db = work.path().join("native.rsdb");
    succeeds(
        tool()
            .args(["database", "build", "--taxonomy"])
            .arg(fixture("taxonomy"))
            .arg("--library")
            .arg(fixture("reference.fna"))
            .arg("--output")
            .arg(&native_db)
            .args([
                "--capacity",
                "1024",
                "--block-size",
                "256",
                "--subblock-size",
                "64",
                "--minimum-bits-for-taxid",
                "28",
                "--threads",
                "4",
            ]),
    );
    let source = std::fs::read(source_db.join("hash.k2d")).unwrap();
    let native = std::fs::read(&native_db).unwrap();
    assert_eq!(&native[native.len() - 4096..], &source[32..]);
    for input in ["reads.fastq", "unknown.fastq"] {
        let source_calls = work.path().join(format!("source-{input}.calls"));
        let source_report = work.path().join(format!("source-{input}.report"));
        succeeds(
            Command::new(oracle.join("kraken2"))
                .arg("--db")
                .arg(&source_db)
                .args(["--threads", "4", "--report-zero-counts"])
                .arg("--output")
                .arg(&source_calls)
                .arg("--report")
                .arg(&source_report)
                .arg(fixture(input)),
        );
        let native_calls = work.path().join(format!("native-{input}.calls"));
        let native_report = work.path().join(format!("native-{input}.report"));
        succeeds(
            tool()
                .args(["classify", "--database"])
                .arg(&native_db)
                .arg("--input")
                .arg(fixture(input))
                .arg("--output")
                .arg(&native_calls)
                .args(["--threads", "4"]),
        );
        succeeds(
            tool()
                .args(["report", "--database"])
                .arg(&native_db)
                .arg("--input")
                .arg(&native_calls)
                .arg("--output")
                .arg(&native_report)
                .arg("--report-zero-counts"),
        );
        assert_eq!(
            std::fs::read(native_calls).unwrap(),
            std::fs::read(source_calls).unwrap()
        );
        assert_eq!(
            std::fs::read(native_report).unwrap(),
            std::fs::read(source_report).unwrap()
        );
    }
}
