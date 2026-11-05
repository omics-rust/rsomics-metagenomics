use std::path::PathBuf;

use clap::{Args, Subcommand};
use rayon::{ThreadPool, ThreadPoolBuilder};
use rsomics_common::{Result, RsomicsError, ThreadArgs};
use rsomics_metagenomics::classifier::{
    self, BuildConfig, ClassifyConfig, InspectConfig, MinimizerOptions, ReportConfig,
};
use serde::Serialize;

#[derive(Debug, Args)]
pub(super) struct DatabaseArgs {
    #[command(subcommand)]
    pub command: DatabaseCommand,
}

#[derive(Debug, Subcommand)]
pub(super) enum DatabaseCommand {
    /// Build a taxonomy-labelled native compact database.
    Build(BuildArgs),
    /// Inspect database parameters and compact-cell taxon counts.
    Inspect(InspectArgs),
}

#[derive(Debug, Args)]
pub(super) struct BuildArgs {
    /// Local NCBI taxdump directory.
    #[arg(long)]
    taxonomy: PathBuf,
    /// Ordered FASTA libraries with explicit taxon IDs or mapped IDs.
    #[arg(long, required = true)]
    library: Vec<PathBuf>,
    /// Optional tab-separated sequence ID/taxon map.
    #[arg(long)]
    mapping: Option<PathBuf>,
    /// New native database file; existing files are not replaced.
    #[arg(short, long)]
    output: PathBuf,
    #[arg(long, default_value_t = 35)]
    kmer_len: usize,
    #[arg(long, default_value_t = 31)]
    minimizer_len: usize,
    #[arg(long, default_value_t = 7)]
    minimizer_spaces: usize,
    /// Explicit compact-cell capacity instead of sampled sizing.
    #[arg(long)]
    capacity: Option<usize>,
    #[arg(long, default_value_t = 16384)]
    block_size: usize,
    #[arg(long, default_value_t = 16384)]
    subblock_size: usize,
    #[arg(long, default_value_t = 0)]
    minimum_bits_for_taxid: u8,
    #[command(flatten)]
    threads: ThreadArgs,
}

#[derive(Debug, Args)]
pub(super) struct InspectArgs {
    #[arg(long)]
    database: PathBuf,
    /// Write the compact-cell taxon report to this new file or -.
    #[arg(short, long)]
    output: Option<PathBuf>,
    #[arg(long)]
    report_zero_counts: bool,
}

#[derive(Debug, Args)]
pub(super) struct ClassifyArgs {
    #[arg(long)]
    database: PathBuf,
    /// Single-end FASTA/FASTQ, with content-detected compression; - reads stdin.
    #[arg(short, long)]
    input: PathBuf,
    /// Five-column classifications; omit or use - for stdout.
    #[arg(short, long)]
    pub output: Option<PathBuf>,
    #[command(flatten)]
    threads: ThreadArgs,
}

#[derive(Debug, Args)]
pub(super) struct ReportArgs {
    #[arg(long)]
    database: PathBuf,
    /// Five-column single-end classifications; - reads stdin.
    #[arg(short, long)]
    input: PathBuf,
    /// Standard six-column report; omit or use - for stdout.
    #[arg(short, long)]
    pub output: Option<PathBuf>,
    #[arg(long)]
    report_zero_counts: bool,
}

#[derive(Serialize)]
#[serde(untagged)]
pub enum Report {
    Build(classifier::DatabaseManifest),
    Inspect(classifier::DatabaseSummary),
    Reads(classifier::ClassifySummary),
}

fn pool(threads: &ThreadArgs) -> Result<ThreadPool> {
    let mut builder = ThreadPoolBuilder::new();
    if let Some(threads) = threads.requested() {
        builder = builder.num_threads(threads.get());
    }
    builder
        .build()
        .map_err(|error| RsomicsError::ConfigError(format!("creating worker pool: {error}")))
}

impl DatabaseArgs {
    pub fn execute(self, json: bool) -> Result<Report> {
        match self.command {
            DatabaseCommand::Build(args) => {
                let pool = pool(&args.threads)?;
                let config = BuildConfig {
                    taxonomy: args.taxonomy,
                    libraries: args.library,
                    mapping: args.mapping,
                    output: args.output,
                    minimizers: MinimizerOptions {
                        k: args.kmer_len,
                        length: args.minimizer_len,
                        spaces: args.minimizer_spaces,
                    },
                    capacity: args.capacity,
                    block_size: args.block_size,
                    subblock_size: args.subblock_size,
                    minimum_taxid_bits: args.minimum_bits_for_taxid,
                };
                classifier::build_database(&config, &pool).map(Report::Build)
            }
            DatabaseCommand::Inspect(args) => {
                if json && args.output.as_deref() == Some(std::path::Path::new("-")) {
                    return Err(RsomicsError::ConfigError(
                        "--json cannot share report stdout".into(),
                    ));
                }
                let output = args.output.or_else(|| (!json).then(|| PathBuf::from("-")));
                classifier::inspect_database(&InspectConfig {
                    database: args.database,
                    output,
                    zero_counts: args.report_zero_counts,
                })
                .map(Report::Inspect)
            }
        }
    }
}

impl ClassifyArgs {
    pub fn execute(self) -> Result<Report> {
        classifier::classify(
            &ClassifyConfig {
                database: self.database,
                input: self.input,
                output: self.output,
            },
            &pool(&self.threads)?,
        )
        .map(Report::Reads)
    }
}

impl ReportArgs {
    pub fn execute(self) -> Result<Report> {
        classifier::report(&ReportConfig {
            database: self.database,
            input: self.input,
            output: self.output,
            zero_counts: self.report_zero_counts,
        })
        .map(Report::Reads)
    }
}
