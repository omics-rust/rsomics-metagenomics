use std::collections::HashMap;
use std::io::Write;

use rsomics_common::Result;
use rsomics_seqio::OwnedRecord;

use super::database::Database;
use super::invalid;
use super::minimizer::Scanner;

const AMBIGUOUS: usize = usize::MAX;

pub(super) fn classify(record: &OwnedRecord, database: &Database) -> Result<(Vec<u8>, bool)> {
    let taxonomy = &database.taxonomy;
    let mut runs: Vec<(usize, usize)> = Vec::new();
    let mut last_key = None;
    let mut last_taxon = 0;
    let mut groups = 0;
    for hit in Scanner::new(
        &record.seq,
        &database.manifest.minimizers,
        0,
        record.seq.len(),
    ) {
        let taxon = if hit.ambiguous {
            AMBIGUOUS
        } else {
            if last_key != Some(hit.key) {
                last_taxon = database.table.get(hit.key);
                last_key = Some(hit.key);
                groups += usize::from(last_taxon != 0);
            }
            last_taxon
        };
        if let Some(last) = runs.last_mut().filter(|last| last.0 == taxon) {
            last.1 += 1;
        } else {
            runs.push((taxon, 1));
        }
    }
    let mut call = 0;
    let mut best = 0;
    if groups >= 2 {
        let mut counts: HashMap<usize, usize> = HashMap::new();
        for &(taxon, count) in &runs {
            if taxon != 0 && taxon != AMBIGUOUS {
                *counts.entry(taxon).or_default() += count;
            }
        }
        for &taxon in counts.keys() {
            let score: usize = counts
                .iter()
                .filter(|&(&ancestor, _)| taxonomy.ancestor(ancestor, taxon))
                .map(|(_, count)| count)
                .sum();
            if score > best {
                call = taxon;
                best = score;
            } else if score == best {
                call = taxonomy.lca(call, taxon);
            }
        }
    }
    let id = record
        .id
        .split(|byte| byte.is_ascii_whitespace())
        .find(|part| !part.is_empty())
        .ok_or_else(|| invalid("empty query ID"))?;
    let mut output = Vec::with_capacity(id.len() + 128);
    output.extend_from_slice(if call == 0 { b"U\t" } else { b"C\t" });
    output.extend_from_slice(id);
    write!(
        output,
        "\t{}\t{}\t",
        taxonomy.nodes[call].external,
        record.seq.len()
    )?;
    if runs.is_empty() {
        output.extend_from_slice(b"0:0");
    }
    for (index, &(taxon, count)) in runs.iter().enumerate() {
        if taxon == AMBIGUOUS {
            write!(output, "A:{count} ")?;
        } else {
            write!(output, "{}:{count}", taxonomy.nodes[taxon].external)?;
            if index + 1 < runs.len() {
                output.push(b' ');
            }
        }
    }
    output.push(b'\n');
    Ok((output, call != 0))
}
