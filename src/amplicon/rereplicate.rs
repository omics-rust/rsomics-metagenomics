use std::io::Write;

use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::OwnedRecord;

use super::{
    CompatibilityProfile, FastaWriter, fasta_record_size, prepare_record, require_abundance,
};

#[derive(Clone, Debug)]
pub struct RereplicateOptions {
    pub profile: CompatibilityProfile,
    pub emit_size_one: bool,
    pub line_width: usize,
    pub maximum_records: u64,
    pub maximum_output_bytes: u64,
}

impl Default for RereplicateOptions {
    fn default() -> Self {
        Self {
            profile: CompatibilityProfile::Strict,
            emit_size_one: false,
            line_width: 0,
            maximum_records: 10_000_000,
            maximum_output_bytes: 1_073_741_824,
        }
    }
}

#[derive(Clone, Debug)]
struct ExpansionRecord {
    id: Vec<u8>,
    sequence: Vec<u8>,
    abundance: u64,
}

#[derive(Clone, Debug)]
pub struct RereplicationPlan {
    records: Vec<ExpansionRecord>,
    output_records: u64,
    output_bytes: u64,
    line_width: usize,
}

impl RereplicationPlan {
    pub fn new<I>(records: I, options: &RereplicateOptions) -> Result<Self>
    where
        I: IntoIterator<Item = OwnedRecord>,
    {
        let mut prepared_records = Vec::new();
        let mut output_records = 0u64;
        let mut output_bytes = 0u64;
        for (input_index, record) in records.into_iter().enumerate() {
            let prepared = prepare_record(record, input_index, options.profile)?;
            let abundance = match options.profile {
                CompatibilityProfile::Strict => require_abundance(&prepared)?,
                CompatibilityProfile::Vsearch => prepared.observed_abundance.unwrap_or(1),
            };
            let mut id = prepared.stripped_id;
            if options.emit_size_one {
                id.extend_from_slice(b";size=1");
            }
            let record_bytes =
                fasta_record_size(id.len(), prepared.sequence.len(), options.line_width)?;
            output_records = output_records.checked_add(abundance).ok_or_else(|| {
                RsomicsError::InvalidInput("rereplication record count overflows u64".into())
            })?;
            output_bytes = record_bytes
                .checked_mul(abundance)
                .and_then(|bytes| output_bytes.checked_add(bytes))
                .ok_or_else(|| {
                    RsomicsError::InvalidInput("rereplication output size overflows u64".into())
                })?;
            prepared_records.push(ExpansionRecord {
                id,
                sequence: prepared.sequence,
                abundance,
            });
        }
        if output_records > options.maximum_records {
            return Err(RsomicsError::ConfigError(format!(
                "rereplication would write {output_records} records, above --maximum-records {}",
                options.maximum_records
            )));
        }
        if output_bytes > options.maximum_output_bytes {
            return Err(RsomicsError::ConfigError(format!(
                "rereplication would write {output_bytes} bytes, above --maximum-output-bytes {}",
                options.maximum_output_bytes
            )));
        }
        Ok(Self {
            records: prepared_records,
            output_records,
            output_bytes,
            line_width: options.line_width,
        })
    }

    #[must_use]
    pub fn input_records(&self) -> usize {
        self.records.len()
    }

    #[must_use]
    pub fn output_records(&self) -> u64 {
        self.output_records
    }

    #[must_use]
    pub fn output_bytes(&self) -> u64 {
        self.output_bytes
    }

    pub fn write<W: Write>(&self, output: W) -> Result<()> {
        let mut writer = FastaWriter::new(output, self.line_width);
        for record in &self.records {
            for _ in 0..record.abundance {
                writer.write_record(&record.id, &record.sequence)?;
            }
        }
        writer.finish().map(drop)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn budget_is_checked_before_writing() {
        let record = OwnedRecord {
            id: b"x;size=3".to_vec(),
            seq: b"ACGT".to_vec(),
            qual: None,
        };
        let options = RereplicateOptions {
            maximum_records: 2,
            ..RereplicateOptions::default()
        };
        assert!(RereplicationPlan::new([record], &options).is_err());
    }

    #[test]
    fn planned_byte_count_matches_written_output() {
        let record = OwnedRecord {
            id: b"x;size=3".to_vec(),
            seq: b"ACGTAC".to_vec(),
            qual: None,
        };
        let options = RereplicateOptions {
            line_width: 4,
            ..RereplicateOptions::default()
        };
        let plan = RereplicationPlan::new([record], &options).unwrap();
        let mut output = Vec::new();
        plan.write(&mut output).unwrap();
        assert_eq!(plan.output_bytes(), output.len() as u64);
    }
}
