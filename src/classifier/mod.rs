mod assign;
mod compact;
mod database;
mod minimizer;
mod report;
mod taxonomy;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use rayon::ThreadPool;
use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError, reject_output_alias};
use rsomics_seqio::{Format, Record};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use compact::CompactTable;
use database::Database;
pub use minimizer::MinimizerOptions;
use minimizer::{Scanner, hash};
use taxonomy::Ncbi;

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}

#[derive(Clone, Debug)]
pub struct BuildConfig {
    pub taxonomy: PathBuf,
    pub libraries: Vec<PathBuf>,
    pub mapping: Option<PathBuf>,
    pub output: PathBuf,
    pub minimizers: MinimizerOptions,
    pub capacity: Option<usize>,
    pub block_size: usize,
    pub subblock_size: usize,
    pub minimum_taxid_bits: u8,
}

pub struct InspectConfig {
    pub database: PathBuf,
    pub output: Option<PathBuf>,
    pub zero_counts: bool,
}

pub struct ClassifyConfig {
    pub database: PathBuf,
    pub input: PathBuf,
    pub output: Option<PathBuf>,
}

pub struct ReportConfig {
    pub database: PathBuf,
    pub input: PathBuf,
    pub output: Option<PathBuf>,
    pub zero_counts: bool,
}

#[derive(Serialize, Deserialize)]
pub struct SourceIdentity {
    pub name: String,
    pub sha256: String,
}

#[derive(Serialize, Deserialize)]
pub struct DatabaseManifest {
    pub format_version: u32,
    pub profile: String,
    pub build_version: String,
    pub minimizers: MinimizerOptions,
    pub capacity: usize,
    pub occupied: usize,
    pub value_bits: u8,
    pub block_size: usize,
    pub subblock_size: usize,
    pub sources: Vec<SourceIdentity>,
}

#[derive(Serialize)]
pub struct DatabaseSummary {
    pub operation: &'static str,
    pub capacity: usize,
    pub occupied: usize,
    pub taxa: usize,
    pub key_bits: u8,
    pub value_bits: u8,
    pub minimizers: MinimizerOptions,
}

#[derive(Serialize)]
pub struct ClassifySummary {
    pub operation: &'static str,
    pub reads: u64,
    pub classified: u64,
    pub unclassified: u64,
}

pub type ReportSummary = ClassifySummary;

fn open_input(path: &Path) -> Result<Box<dyn Read>> {
    if path == Path::new("-") {
        Ok(Box::new(io::stdin()))
    } else {
        Ok(Box::new(File::open(path)?))
    }
}

fn staged_output(path: &Path) -> Result<tempfile::NamedTempFile> {
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let mut builder = tempfile::Builder::new();
    builder.prefix(".rsomics-");
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        builder.permissions(std::fs::Permissions::from_mode(0o666));
    }
    Ok(builder.tempfile_in(parent)?)
}

fn write_output<T>(
    path: Option<&Path>,
    inputs: &[&Path],
    operation: impl FnOnce(&mut dyn Write) -> Result<T>,
) -> Result<T> {
    if let Some(path) = path.filter(|path| *path != Path::new("-")) {
        reject_output_alias(path, inputs.iter().copied())?;
        if path.try_exists()? {
            return Err(invalid(format!("output {} already exists", path.display())));
        }
        let mut temporary = staged_output(path)?;
        let result = {
            let mut writer = BufWriter::with_capacity(1024 * 1024, temporary.as_file_mut());
            let result = operation(&mut writer)?;
            writer.flush()?;
            result
        };
        temporary
            .persist_noclobber(path)
            .map_err(|error| RsomicsError::Io(error.error))?;
        Ok(result)
    } else {
        rsomics_common::write_output(None, operation)
    }
}

fn source_identity(path: &Path) -> Result<SourceIdentity> {
    let mut file = File::open(path)?;
    let mut digest = Sha256::new();
    let mut buffer = [0u8; 64 * 1024];
    loop {
        let bytes = file.read(&mut buffer)?;
        if bytes == 0 {
            break;
        }
        digest.update(&buffer[..bytes]);
    }
    Ok(SourceIdentity {
        name: path.display().to_string(),
        sha256: format!("{:x}", digest.finalize()),
    })
}

fn id_taxon(id: &[u8], mapping: &HashMap<String, u64>) -> Result<u64> {
    let id = id
        .split(|byte| byte.is_ascii_whitespace())
        .find(|part| !part.is_empty())
        .ok_or_else(|| invalid("empty reference ID"))?;
    let id = std::str::from_utf8(id).map_err(|_| invalid("reference ID is not UTF-8"))?;
    let mut parts = id.split('|');
    while let Some(part) = parts.next() {
        if part == "kraken:taxid" {
            return parts
                .next()
                .ok_or_else(|| invalid("missing explicit taxon ID"))?
                .parse()
                .map_err(|_| invalid("invalid explicit taxon ID"));
        }
    }
    mapping
        .get(id)
        .copied()
        .ok_or_else(|| invalid(format!("reference {id} has no taxon assignment")))
}

fn load_mapping(path: Option<&Path>) -> Result<HashMap<String, u64>> {
    let mut mapping = HashMap::new();
    if let Some(path) = path {
        for line in BufReader::new(File::open(path)?).lines() {
            let line = line?;
            let (id, taxon) = line
                .split_once('\t')
                .ok_or_else(|| invalid("mapping must contain ID and taxon separated by one tab"))?;
            let taxon = taxon
                .parse::<u64>()
                .map_err(|_| invalid("invalid mapped taxon ID"))?;
            if id.is_empty()
                || id.bytes().any(|b| b.is_ascii_whitespace())
                || mapping.insert(id.to_owned(), taxon).is_some()
            {
                return Err(invalid(
                    "empty, whitespace-containing or duplicate mapped ID",
                ));
            }
        }
    }
    Ok(mapping)
}

fn visit_references(
    config: &BuildConfig,
    mut visitor: impl FnMut(Record<'_>) -> Result<()>,
) -> Result<()> {
    for path in &config.libraries {
        let mut reader = rsomics_seqio::open_path(path)?;
        if reader.format() != Format::Fasta {
            return Err(invalid("database libraries must be FASTA"));
        }
        while let Some(record) = reader.read_record()? {
            visitor(record)?;
        }
    }
    Ok(())
}

fn block_keys(sequence: &[u8], start: usize, config: &BuildConfig) -> Vec<u64> {
    let finish = sequence.len().min(
        start
            .saturating_add(config.block_size)
            .saturating_add(config.minimizers.k - 1),
    );
    let mut keys = Vec::new();
    for part in (start..finish).step_by(config.subblock_size) {
        let end = finish.min(
            part.saturating_add(config.subblock_size)
                .saturating_add(config.minimizers.k - 1),
        );
        keys.extend(
            Scanner::new(sequence, &config.minimizers, part, end)
                .filter(|hit| !hit.ambiguous)
                .map(|hit| hit.key),
        );
    }
    keys.sort_unstable_by_key(|&key| (hash(key) & 255, key));
    keys.dedup();
    keys
}

pub fn build_database(config: &BuildConfig, pool: &ThreadPool) -> Result<DatabaseManifest> {
    config.minimizers.validate()?;
    if config.output == Path::new("-")
        || config.libraries.iter().any(|path| path == Path::new("-"))
        || config.libraries.is_empty()
        || config.block_size == 0
        || config.subblock_size == 0
        || config.subblock_size > config.block_size
        || config.minimum_taxid_bits > 31
    {
        return Err(invalid(
            "invalid database libraries, block/subblock sizes or taxon bit width",
        ));
    }
    let mut inputs: Vec<&Path> = config.libraries.iter().map(PathBuf::as_path).collect();
    let taxonomy_files: Vec<_> = ["nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"]
        .iter()
        .map(|name| config.taxonomy.join(name))
        .collect();
    for path in &taxonomy_files {
        if path.try_exists()? {
            inputs.push(path);
        }
    }
    if let Some(mapping) = &config.mapping {
        inputs.push(mapping);
    }
    reject_output_alias(&config.output, inputs.iter().copied())?;
    if config.output.try_exists()? {
        return Err(invalid("database output already exists"));
    }
    let sources = inputs
        .iter()
        .map(|path| source_identity(path))
        .collect::<Result<Vec<_>>>()?;
    let ncbi = Ncbi::load(&config.taxonomy)?;
    let mapping = load_mapping(config.mapping.as_deref())?;
    let mut taxa = HashSet::new();
    let mut assignments = HashMap::new();
    let mut samples = HashSet::new();
    let mut records = 0u64;
    visit_references(config, |record| {
        let original = id_taxon(record.id, &mapping)?;
        let resolved = ncbi.resolve(original)?;
        assignments.insert(original, resolved);
        taxa.insert(resolved);
        for hit in Scanner::new(record.seq, &config.minimizers, 0, record.seq.len()) {
            if !hit.ambiguous && (hash(hit.key) & 1023) < 4 {
                samples.insert(hit.key);
            }
        }
        records += 1;
        Ok(())
    })?;
    if records == 0 {
        return Err(invalid("database libraries contain no sequences"));
    }
    let taxonomy = ncbi.prune(&taxa)?;
    drop(ncbi);
    drop(taxa);
    let needed = (usize::BITS - (taxonomy.nodes.len() - 1).leading_zeros()) as u8;
    let bits = needed.max(config.minimum_taxid_bits);
    if bits > 31 {
        return Err(invalid("taxonomy cannot fit compact values"));
    }
    let estimate = samples
        .len()
        .checked_add(1)
        .and_then(|n| n.checked_mul(256))
        .and_then(|n| n.checked_add(8192))
        .ok_or_else(|| invalid("capacity estimate overflow"))?;
    let capacity = config.capacity.unwrap_or((estimate as f64 / 0.7) as usize);
    drop(samples);
    capacity
        .checked_mul(4)
        .ok_or_else(|| invalid("database capacity overflow"))?;
    let mut table = CompactTable::new(capacity, bits)?;
    visit_references(config, |record| {
        let original = id_taxon(record.id, &mapping)?;
        let resolved = assignments
            .get(&original)
            .copied()
            .ok_or_else(|| invalid("reference taxon changed during build"))?;
        let taxon = taxonomy.internal(resolved)?;
        let mut starts = (0..record.seq.len()).step_by(config.block_size);
        loop {
            let batch: Vec<_> = starts.by_ref().take(32).collect();
            if batch.is_empty() {
                break;
            }
            let blocks: Vec<_> = pool.install(|| {
                batch
                    .par_iter()
                    .map(|&start| block_keys(record.seq, start, config))
                    .collect()
            });
            for keys in blocks {
                for key in keys {
                    table.insert(key, taxon, &taxonomy)?;
                }
            }
        }
        Ok(())
    })?;
    for (path, identity) in inputs.iter().zip(&sources) {
        if source_identity(path)?.sha256 != identity.sha256 {
            return Err(invalid("database source changed during build"));
        }
    }
    let manifest = DatabaseManifest {
        format_version: 1,
        profile: "kraken2-2.17.1-linear-v1".into(),
        build_version: env!("CARGO_PKG_VERSION").into(),
        minimizers: config.minimizers.clone(),
        capacity,
        occupied: table.occupied,
        value_bits: bits,
        block_size: config.block_size,
        subblock_size: config.subblock_size,
        sources,
    };
    database::save(&config.output, &manifest, &taxonomy, &table)?;
    Ok(manifest)
}

pub fn inspect_database(config: &InspectConfig) -> Result<DatabaseSummary> {
    let database = Database::load(&config.database)?;
    let summary = DatabaseSummary {
        operation: "database-inspect",
        capacity: database.table.cells.len(),
        occupied: database.table.occupied,
        taxa: database.taxonomy.nodes.len() - 1,
        key_bits: 32 - database.table.value_bits,
        value_bits: database.table.value_bits,
        minimizers: database.manifest.minimizers.clone(),
    };
    if let Some(output) = &config.output {
        let mut counts = report::Counts::new(database.taxonomy.nodes.len());
        for &cell in &database.table.cells {
            let taxon = database.table.value(cell);
            if taxon != 0 {
                counts.add(taxon)?;
            }
        }
        write_output(Some(output), &[&config.database], |writer| {
            counts.write(writer, &database.taxonomy, config.zero_counts)
        })?;
    }
    Ok(summary)
}

pub fn classify(config: &ClassifyConfig, pool: &ThreadPool) -> Result<ClassifySummary> {
    let database = Database::load(&config.database)?;
    let mut reader = rsomics_seqio::open_reader(open_input(&config.input)?)?;
    write_output(
        config.output.as_deref(),
        &[&config.database, &config.input],
        |writer| {
            let mut summary = ClassifySummary {
                operation: "classify",
                reads: 0,
                classified: 0,
                unclassified: 0,
            };
            loop {
                let mut batch = Vec::new();
                let mut bytes = 0;
                while bytes < 1024 * 1024 {
                    let Some(record) = reader.read_record()? else {
                        break;
                    };
                    bytes += record.seq.len() + record.id.len();
                    batch.push(record.to_owned());
                }
                if batch.is_empty() {
                    break;
                }
                let results: Vec<_> = pool.install(|| {
                    batch
                        .par_iter()
                        .map(|record| assign::classify(record, &database))
                        .collect()
                });
                for result in results {
                    let (bytes, classified) = result?;
                    writer.write_all(&bytes)?;
                    summary.reads += 1;
                    summary.classified += u64::from(classified);
                    summary.unclassified += u64::from(!classified);
                }
            }
            Ok(summary)
        },
    )
}

pub fn report(config: &ReportConfig) -> Result<ReportSummary> {
    let database = Database::load(&config.database)?;
    let mut counts = report::Counts::new(database.taxonomy.nodes.len());
    report::read_calls(
        BufReader::new(open_input(&config.input)?),
        &database.taxonomy,
        &mut counts,
    )?;
    write_output(
        config.output.as_deref(),
        &[&config.database, &config.input],
        |writer| counts.write(writer, &database.taxonomy, config.zero_counts),
    )?;
    Ok(ClassifySummary {
        operation: "report",
        reads: counts.total,
        classified: counts.total - counts.direct[0],
        unclassified: counts.direct[0],
    })
}
