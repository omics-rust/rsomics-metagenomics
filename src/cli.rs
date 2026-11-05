use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

mod classifier;

use clap::{Args, Parser, Subcommand, ValueEnum};
use rsomics_common::{
    Context, OutputArgs, Result, RsomicsError, ToolMeta, reject_output_alias, write_output,
};
use rsomics_metagenomics::amplicon::{
    Amplicon, CompatibilityProfile, DereplicateMode, DereplicateOptions, Dereplication,
    FastaWriter, RereplicateOptions, RereplicationPlan, SortOptions, dereplicate_fasta, read_fasta,
    sort_by_abundance,
};
use rsomics_seqio::OwnedRecord;
use serde::Serialize;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-metagenomics",
    version,
    about = "Run amplicon and taxonomic classification workflows",
    arg_required_else_help = true
)]
pub struct Cli {
    #[command(subcommand)]
    command: Command,

    #[command(flatten)]
    pub output: OutputArgs,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Collapse identical or prefix-related amplicons.
    Dereplicate(DereplicateArgs),
    /// Sort amplicons by decreasing size abundance.
    SortAbundance(SortArgs),
    /// Expand abundance-labelled amplicons into individual records.
    Rereplicate(RereplicateArgs),
    /// Build or inspect a native taxonomic database.
    Database(classifier::DatabaseArgs),
    /// Classify single-end reads against a native database.
    Classify(classifier::ClassifyArgs),
    /// Summarize classifications as a standard taxonomic report.
    Report(classifier::ReportArgs),
}

impl Command {
    fn output(&self) -> Option<&Path> {
        match self {
            Self::Dereplicate(args) => args.io.output.as_deref(),
            Self::SortAbundance(args) => args.io.output.as_deref(),
            Self::Rereplicate(args) => args.io.output.as_deref(),
            Self::Database(_) => None,
            Self::Classify(args) => args.output.as_deref(),
            Self::Report(args) => args.output.as_deref(),
        }
    }
}

#[derive(Debug, Args)]
struct InputOutput {
    /// Input FASTA path; use - for standard input.
    #[arg(short, long)]
    input: PathBuf,

    /// Transactional FASTA output; omit or use - for standard output.
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
enum ProfileArg {
    #[default]
    Strict,
    Vsearch,
}

impl From<ProfileArg> for CompatibilityProfile {
    fn from(value: ProfileArg) -> Self {
        match value {
            ProfileArg::Strict => Self::Strict,
            ProfileArg::Vsearch => Self::Vsearch,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
enum DereplicateModeArg {
    #[default]
    FullLength,
    Prefix,
}

impl From<DereplicateModeArg> for DereplicateMode {
    fn from(value: DereplicateModeArg) -> Self {
        match value {
            DereplicateModeArg::FullLength => Self::FullLength,
            DereplicateModeArg::Prefix => Self::Prefix,
        }
    }
}

#[derive(Debug, Args)]
struct DereplicateArgs {
    #[command(flatten)]
    io: InputOutput,

    /// Sequence relationship used to form clusters.
    #[arg(long, value_enum, default_value_t)]
    mode: DereplicateModeArg,

    /// Header and edge-case behavior profile.
    #[arg(long, value_enum, default_value_t)]
    profile: ProfileArg,

    /// Sum input size attributes instead of counting every record once.
    #[arg(long)]
    input_abundance: bool,

    /// Minimum accepted sequence length.
    #[arg(long, default_value_t = 32)]
    minimum_length: usize,

    /// Maximum accepted sequence length.
    #[arg(long, default_value_t = 50_000)]
    maximum_length: usize,

    /// Minimum emitted cluster abundance.
    #[arg(long, default_value_t = 1)]
    minimum_abundance: u64,

    /// Maximum emitted cluster abundance.
    #[arg(long)]
    maximum_abundance: Option<u64>,

    /// Emit only the first N clusters after filtering and ordering.
    #[arg(long)]
    top: Option<usize>,

    /// Sequence line width; defaults to 0 in strict and 80 in vsearch.
    #[arg(long)]
    line_width: Option<usize>,
}

#[derive(Debug, Args)]
struct SortArgs {
    #[command(flatten)]
    io: InputOutput,

    /// Header and edge-case behavior profile.
    #[arg(long, value_enum, default_value_t)]
    profile: ProfileArg,

    /// Minimum accepted sequence length.
    #[arg(long, default_value_t = 1)]
    minimum_length: usize,

    /// Maximum accepted sequence length.
    #[arg(long, default_value_t = 50_000)]
    maximum_length: usize,

    /// Minimum accepted abundance.
    #[arg(long, default_value_t = 1)]
    minimum_abundance: u64,

    /// Maximum accepted abundance.
    #[arg(long)]
    maximum_abundance: Option<u64>,

    /// Emit only the first N records after filtering and ordering.
    #[arg(long)]
    top: Option<usize>,

    /// Sequence line width; defaults to 0 in strict and 80 in vsearch.
    #[arg(long)]
    line_width: Option<usize>,
}

#[derive(Debug, Args)]
struct RereplicateArgs {
    #[command(flatten)]
    io: InputOutput,

    /// Header and edge-case behavior profile.
    #[arg(long, value_enum, default_value_t)]
    profile: ProfileArg,

    /// Add size=1 to every expanded record.
    #[arg(long)]
    emit_size_one: bool,

    /// Sequence line width; defaults to 0 in strict and 80 in vsearch.
    #[arg(long)]
    line_width: Option<usize>,

    /// Refuse expansions above this record count.
    #[arg(long, default_value_t = 10_000_000)]
    maximum_records: u64,

    /// Refuse expansions above this exact FASTA byte count.
    #[arg(long, default_value_t = 1_073_741_824)]
    maximum_output_bytes: u64,
}

#[derive(Debug, Serialize)]
pub struct AmpliconReport {
    operation: &'static str,
    input_records: usize,
    output_records: u64,
    output_bytes: Option<u64>,
    discarded_records: usize,
    output: String,
}

#[derive(Serialize)]
#[serde(untagged)]
pub enum ExecutionReport {
    Amplicon(AmpliconReport),
    Classifier(classifier::Report),
}

impl Cli {
    pub fn execute(self) -> Result<ExecutionReport> {
        if self.output.json
            && !matches!(self.command, Command::Database(_))
            && is_stdout(self.command.output())
        {
            return Err(RsomicsError::ConfigError(
                "--json requires a named data output".into(),
            ));
        }
        match self.command {
            Command::Dereplicate(args) => execute_dereplicate(args).map(ExecutionReport::Amplicon),
            Command::SortAbundance(args) => execute_sort(args).map(ExecutionReport::Amplicon),
            Command::Rereplicate(args) => execute_rereplicate(args).map(ExecutionReport::Amplicon),
            Command::Database(args) => args
                .execute(self.output.json)
                .map(ExecutionReport::Classifier),
            Command::Classify(args) => args.execute().map(ExecutionReport::Classifier),
            Command::Report(args) => args.execute().map(ExecutionReport::Classifier),
        }
    }
}

fn execute_dereplicate(args: DereplicateArgs) -> Result<AmpliconReport> {
    reject_alias(&args.io)?;
    let profile = args.profile.into();
    let options = DereplicateOptions {
        mode: args.mode.into(),
        profile,
        input_abundance: args.input_abundance,
        minimum_length: args.minimum_length,
        maximum_length: args.maximum_length,
        minimum_abundance: args.minimum_abundance,
        maximum_abundance: args.maximum_abundance.unwrap_or(u64::MAX),
        top: args.top,
    };
    let result = dereplicate_input(&args.io.input, &options)?;
    let width = line_width(args.line_width, profile);
    write_amplicons(result.records(), args.io.output.as_deref(), width)?;
    Ok(AmpliconReport {
        operation: "dereplicate",
        input_records: result.input_records(),
        output_records: result.records().len() as u64,
        output_bytes: None,
        discarded_records: result.discarded_length()
            + result
                .clusters_before_filter()
                .saturating_sub(result.records().len()),
        output: output_label(args.io.output.as_deref()),
    })
}

fn execute_sort(args: SortArgs) -> Result<AmpliconReport> {
    reject_alias(&args.io)?;
    let profile = args.profile.into();
    let records = read_input(&args.io.input, profile)?;
    let result = sort_by_abundance(
        records,
        &SortOptions {
            profile,
            minimum_length: args.minimum_length,
            maximum_length: args.maximum_length,
            minimum_abundance: args.minimum_abundance,
            maximum_abundance: args.maximum_abundance.unwrap_or(u64::MAX),
            top: args.top,
        },
    )?;
    let width = line_width(args.line_width, profile);
    write_amplicons(result.records(), args.io.output.as_deref(), width)?;
    Ok(AmpliconReport {
        operation: "sort-abundance",
        input_records: result.input_records(),
        output_records: result.records().len() as u64,
        output_bytes: None,
        discarded_records: result
            .input_records()
            .saturating_sub(result.records().len()),
        output: output_label(args.io.output.as_deref()),
    })
}

fn execute_rereplicate(args: RereplicateArgs) -> Result<AmpliconReport> {
    reject_alias(&args.io)?;
    let profile = args.profile.into();
    let records = read_input(&args.io.input, profile)?;
    let plan = RereplicationPlan::new(
        records,
        &RereplicateOptions {
            profile,
            emit_size_one: args.emit_size_one,
            line_width: line_width(args.line_width, profile),
            maximum_records: args.maximum_records,
            maximum_output_bytes: args.maximum_output_bytes,
        },
    )?;
    write_output(args.io.output.as_deref(), |writer| plan.write(writer))?;
    Ok(AmpliconReport {
        operation: "rereplicate",
        input_records: plan.input_records(),
        output_records: plan.output_records(),
        output_bytes: Some(plan.output_bytes()),
        discarded_records: 0,
        output: output_label(args.io.output.as_deref()),
    })
}

fn read_input(path: &Path, profile: CompatibilityProfile) -> Result<Vec<OwnedRecord>> {
    let allow_empty = profile == CompatibilityProfile::Vsearch;
    if path == Path::new("-") {
        read_fasta(io::stdin().lock(), allow_empty)
    } else {
        let input = File::open(path).rs_with_context(|| format!("opening {}", path.display()))?;
        read_fasta(input, allow_empty).rs_with_context(|| format!("reading {}", path.display()))
    }
}

fn dereplicate_input(path: &Path, options: &DereplicateOptions) -> Result<Dereplication> {
    if path == Path::new("-") {
        dereplicate_fasta(io::stdin().lock(), options)
    } else {
        let input = File::open(path).rs_with_context(|| format!("opening {}", path.display()))?;
        dereplicate_fasta(input, options).rs_with_context(|| format!("reading {}", path.display()))
    }
}

fn write_amplicons(records: &[Amplicon], output: Option<&Path>, line_width: usize) -> Result<()> {
    write_output(output, |writer| {
        let mut fasta = FastaWriter::new(writer, line_width);
        for record in records {
            let mut id = record.id().to_vec();
            write!(&mut id, ";size={}", record.abundance())?;
            fasta.write_record(&id, record.sequence())?;
        }
        fasta.finish().map(drop)
    })
}

fn reject_alias(io: &InputOutput) -> Result<()> {
    if let Some(output) = io.output.as_deref() {
        reject_output_alias(output, [io.input.as_path()])?;
    }
    Ok(())
}

fn line_width(requested: Option<usize>, profile: CompatibilityProfile) -> usize {
    requested.unwrap_or(if profile == CompatibilityProfile::Vsearch {
        80
    } else {
        0
    })
}

fn is_stdout(output: Option<&Path>) -> bool {
    output.is_none() || output == Some(Path::new("-"))
}

fn output_label(output: Option<&Path>) -> String {
    output.map_or_else(|| "stdout".into(), |path| path.display().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn command_tree_is_valid() {
        Cli::command().debug_assert();
    }
}
