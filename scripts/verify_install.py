import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess


def identity(path):
    with path.open("rb") as stream:
        return {"bytes": path.stat().st_size, "sha256": hashlib.file_digest(stream, "sha256").hexdigest()}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def classify(run, binary, kraken, work, fixture):
    taxonomy = fixture / "taxonomy"
    reference = fixture / "reference.fna"
    reads = fixture / "reads.fastq"
    identifiers = [line[1:] for line in reference.read_bytes().splitlines() if line.startswith(b">")]
    mapping = work / "mapping.tsv"
    mapping.write_bytes(b"".join(name + b"\t" + name.rsplit(b"|", 1)[1] + b"\n" for name in identifiers))
    source_db = work / "source-db"
    source_db.mkdir()
    run(kraken / "build_db", "-k", "35", "-l", "31", "-S", "1111111111111111101010101010101",
        "-c", "1024", "-p", "4", "-B", "256", "-b", "64", "-r", "28", "-n", taxonomy,
        "-m", mapping, "-H", source_db / "hash.k2d", "-t", source_db / "taxo.k2d",
        "-o", source_db / "opts.k2d", stdin=reference.read_bytes())
    database = work / "native.rsdb"
    run(binary, "database", "build", "--taxonomy", taxonomy, "--library", reference,
        "--output", database, "--capacity", "1024", "--block-size", "256", "--subblock-size", "64",
        "--minimum-bits-for-taxid", "28", "--threads", "4")
    result = json.loads(run(binary, "--json", "database", "inspect", "--database", database).stdout)["result"]
    require(result["capacity"] == 1024 and result["occupied"] == 662, "unexpected installed database geometry")
    run(kraken / "kraken2", "--db", source_db, "--threads", "4", "--report-zero-counts",
        "--output", "source-calls.tsv", "--report", "source-report.tsv", reads)
    run(binary, "classify", "--database", database, "--input", reads, "--output", "calls.tsv", "--threads", "4")
    run(binary, "report", "--database", database, "--input", "calls.tsv", "--output", "report.tsv",
        "--report-zero-counts")
    for name in ("calls.tsv", "report.tsv"):
        expected = (work / f"source-{name}").read_bytes()
        require(expected and (work / name).read_bytes() == expected, f"installed {name} differs from Kraken")


def amplicon(run, binary, vsearch, work):
    native_input = source_input = "amplicon.fasta"
    (work / native_input).write_bytes(
        b">first;size=2\n" + b"ACGT" * 10 + b"\n>duplicate;size=4\n" + b"ACGT" * 10
        + b"\n>second;size=3\n" + b"TTGC" * 10 + b"\n"
    )
    for command, oracle_flag, native_options, oracle_options, records in (
        ("dereplicate", "--derep_fulllength", ("--mode", "full-length", "--input-abundance"), ("--sizein",), 2),
        ("sort-abundance", "--sortbysize", ("--minimum-abundance", "4"), ("--minsize", "4"), 1),
        ("rereplicate", "--rereplicate", ("--emit-size-one",), (), 6),
    ):
        native_output, source_output = f"{command}.fasta", f"source-{command}.fasta"
        run(binary, command, "--profile", "vsearch", "--input", native_input, "--output", native_output, *native_options)
        run(vsearch, oracle_flag, source_input, "--output", source_output, "--sizeout", "--quiet", *oracle_options)
        expected = (work / source_output).read_bytes()
        require(sum(line.startswith(b">") for line in expected.splitlines()) == records,
                f"unexpected VSEARCH record count for {command}")
        require((work / native_output).read_bytes() == expected, f"installed {command} output differs from VSEARCH")
        native_input, source_input = native_output, source_output


def verify(binary, kraken, vsearch, delivery):
    source = json.loads((delivery / "SOURCE.json").read_bytes())
    work = delivery / "lifecycle"
    work.mkdir()
    calls = []

    def save(name, value):
        (work / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")

    def run(program, *args, stdin=None):
        argv = list(map(str, (program, *args)))
        result = subprocess.run(argv, cwd=work, input=stdin, capture_output=True, timeout=120,
                                env=dict(os.environ, NO_COLOR="1"))
        calls.append({"argv": argv, "returncode": result.returncode})
        (work / f"{len(calls)}.stdout").write_bytes(result.stdout)
        (work / f"{len(calls)}.stderr").write_bytes(result.stderr)
        save("commands.json", calls)
        require(result.returncode == 0, f"exit {result.returncode}: {argv}; {result.stderr.decode(errors='replace')}")
        return result

    require(source["product"] == "rsomics-metagenomics", "incorrect source package")
    version = f"rsomics-metagenomics {source['version']}".encode()
    require(run(binary, "--version").stdout.strip() == version, "incorrect installed version")
    require(run(kraken / "kraken2", "--version").stdout.startswith(b"Kraken version 2.17.1\n"),
            "expected Kraken 2.17.1")
    upstream = run(vsearch, "--version")
    require(b"vsearch v2.31.0" in upstream.stdout + upstream.stderr, "expected VSEARCH 2.31.0")
    classify(run, binary, kraken, work, Path(__file__).resolve().parents[1] / "tests/data/classification")
    amplicon(run, binary, vsearch, work)
    save("VERIFIED.json", {
        "binary": identity(binary), "source_head": source["head_sha"], "version": source["version"],
        "native_calls": sum(row["argv"][0] == str(binary) for row in calls),
        "kraken": {"version": "2.17.1", "launcher": identity(kraken / "kraken2"),
                   "build_db": identity(kraken / "build_db"), "classify": identity(kraken / "classify")},
        "vsearch": {"version": "2.31.0", "binary": identity(vsearch)},
        "classification_and_amplicon_outputs_match": True,
    })


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for option in ("binary", "kraken-dir", "vsearch", "delivery"):
        parser.add_argument("--" + option, type=Path, required=True)
    options = parser.parse_args()
    verify(*(p.resolve(strict=True) for p in (options.binary, options.kraken_dir, options.vsearch, options.delivery)))
