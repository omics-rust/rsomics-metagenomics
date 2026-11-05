use std::io::{BufRead, Write};

use rsomics_common::Result;

use super::{invalid, taxonomy::Taxonomy};

pub(super) struct Counts {
    pub direct: Vec<u64>,
    pub total: u64,
    observed: Vec<bool>,
}

impl Counts {
    pub fn new(nodes: usize) -> Self {
        Self {
            direct: vec![0; nodes],
            total: 0,
            observed: vec![false; nodes],
        }
    }

    pub fn add(&mut self, taxon: usize) -> Result<()> {
        self.direct[taxon] = self.direct[taxon]
            .checked_add(1)
            .ok_or_else(|| invalid("report direct count overflow"))?;
        self.total = self
            .total
            .checked_add(1)
            .ok_or_else(|| invalid("report total overflow"))?;
        self.observed[taxon] = true;
        Ok(())
    }

    pub fn write(
        &self,
        output: &mut dyn Write,
        taxonomy: &Taxonomy,
        zero_counts: bool,
    ) -> Result<()> {
        if self.total == 0 {
            return Err(invalid(
                "cannot calculate report percentages for zero input sequences",
            ));
        }
        let mut clade = self.direct.clone();
        let mut present = self.observed.clone();
        for index in (1..clade.len()).rev() {
            let parent = taxonomy.nodes[index].parent;
            if parent != 0 {
                clade[parent] = clade[parent]
                    .checked_add(clade[index])
                    .ok_or_else(|| invalid("report clade count overflow"))?;
                present[parent] |= present[index];
            }
        }
        if self.direct[0] != 0 || zero_counts {
            writeln!(
                output,
                "{:6.2}\t{}\t{}\tU\t0\tunclassified",
                100.0 * self.direct[0] as f64 / self.total as f64,
                self.direct[0],
                self.direct[0]
            )?;
        }
        let mut stack = vec![(1usize, 'R', -1i32, 0usize)];
        while let Some((id, inherited, rank_depth, depth)) = stack.pop() {
            if clade[id] == 0 && !zero_counts {
                continue;
            }
            let node = &taxonomy.nodes[id];
            let standard = match node.rank.as_str() {
                "superkingdom" | "domain" => Some('D'),
                "kingdom" => Some('K'),
                "phylum" => Some('P'),
                "class" => Some('C'),
                "order" => Some('O'),
                "family" => Some('F'),
                "genus" => Some('G'),
                "species" => Some('S'),
                _ => None,
            };
            let (code, rank_depth) = standard.map_or((inherited, rank_depth + 1), |code| (code, 0));
            let rank = if rank_depth == 0 {
                code.to_string()
            } else {
                format!("{code}{rank_depth}")
            };
            writeln!(
                output,
                "{:6.2}\t{}\t{}\t{}\t{}\t{}{}",
                100.0 * clade[id] as f64 / self.total as f64,
                clade[id],
                self.direct[id],
                rank,
                node.external,
                "  ".repeat(depth),
                node.name
            )?;
            let mut children = taxonomy.children[id].clone();
            children.sort_by(|&first, &second| {
                present[second]
                    .cmp(&present[first])
                    .then_with(|| clade[second].cmp(&clade[first]))
                    .then_with(|| {
                        taxonomy.nodes[first]
                            .external
                            .cmp(&taxonomy.nodes[second].external)
                    })
            });
            for child in children.into_iter().rev() {
                stack.push((child, code, rank_depth, depth + 1));
            }
        }
        Ok(())
    }
}

pub(super) fn read_calls(
    input: impl BufRead,
    taxonomy: &Taxonomy,
    counts: &mut Counts,
) -> Result<()> {
    for (index, line) in input.lines().enumerate() {
        let line = line?;
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 5
            || fields[1].is_empty()
            || fields[1].bytes().any(|b| b.is_ascii_whitespace())
        {
            return Err(invalid(format!(
                "classification row {} must have five single-end fields",
                index + 1
            )));
        }
        let external = fields[2]
            .parse::<u64>()
            .map_err(|_| invalid("invalid classification taxon"))?;
        if !matches!((fields[0], external), ("U", 0) | ("C", 1..)) {
            return Err(invalid("classification status/taxon mismatch"));
        }
        fields[3]
            .parse::<usize>()
            .map_err(|_| invalid("invalid single-end classification length"))?;
        let taxon = taxonomy.internal(external)?;
        let mut groups = 0usize;
        for token in fields[4].split_ascii_whitespace() {
            let (code, count) = token
                .split_once(':')
                .ok_or_else(|| invalid("invalid classification hit list"))?;
            let count = count
                .parse::<usize>()
                .map_err(|_| invalid("invalid hit count"))?;
            if code != "A" {
                let code = code
                    .parse::<u64>()
                    .map_err(|_| invalid("invalid hit taxon"))?;
                let hit = taxonomy.internal(code)?;
                if hit != 0 && count > 0 {
                    counts.observed[hit] = true;
                }
            }
            groups += 1;
        }
        if groups == 0 {
            return Err(invalid("empty classification hit list"));
        }
        counts.add(taxon)?;
    }
    Ok(())
}
