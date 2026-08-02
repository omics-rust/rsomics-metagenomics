mod abundance;
mod dereplicate;
mod io;
mod rereplicate;
mod sort;

use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::OwnedRecord;

pub use dereplicate::{
    DereplicateMode, DereplicateOptions, Dereplication, dereplicate, dereplicate_fasta,
};
pub use io::{FastaWriter, fasta_record_size, read_fasta, visit_fasta};
pub use rereplicate::{RereplicateOptions, RereplicationPlan};
pub use sort::{AbundanceSort, SortOptions, sort_by_abundance};

use abundance::parse_header;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum CompatibilityProfile {
    #[default]
    Strict,
    Vsearch,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Amplicon {
    id: Vec<u8>,
    sequence: Vec<u8>,
    abundance: u64,
    sort_key: Vec<u8>,
    input_index: usize,
    hash_next: Option<usize>,
}

impl Amplicon {
    #[must_use]
    pub fn id(&self) -> &[u8] {
        &self.id
    }

    #[must_use]
    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }

    #[must_use]
    pub fn abundance(&self) -> u64 {
        self.abundance
    }
}

#[derive(Clone, Debug)]
struct PreparedRecord {
    raw_id: Vec<u8>,
    stripped_id: Vec<u8>,
    sequence: Vec<u8>,
    observed_abundance: Option<u64>,
    input_index: usize,
}

fn prepare_record(
    record: OwnedRecord,
    input_index: usize,
    profile: CompatibilityProfile,
) -> Result<PreparedRecord> {
    if record.qual.is_some() {
        return Err(RsomicsError::InvalidInput(
            "amplicon operations require FASTA input".into(),
        ));
    }
    if record.id.is_empty()
        || record
            .id
            .iter()
            .any(|byte| *byte != b'\t' && !(b' '..=b'~').contains(byte))
    {
        return Err(RsomicsError::InvalidInput(
            "FASTA identifier is empty or contains an invalid byte".into(),
        ));
    }
    if record.seq.is_empty() || record.seq.iter().any(|byte| !(b'!'..=b'~').contains(byte)) {
        return Err(RsomicsError::InvalidInput(
            "FASTA sequence is empty or contains an invalid byte".into(),
        ));
    }
    let parsed = parse_header(record.id, profile)?;
    Ok(PreparedRecord {
        stripped_id: parsed.stripped(),
        raw_id: parsed.raw,
        sequence: record.seq,
        observed_abundance: parsed.abundance,
        input_index,
    })
}

fn require_abundance(record: &PreparedRecord) -> Result<u64> {
    record.observed_abundance.ok_or_else(|| {
        RsomicsError::InvalidInput(format!(
            "record {:?} has no size abundance",
            String::from_utf8_lossy(&record.raw_id)
        ))
    })
}

fn checked_sum(left: u64, right: u64, id: &[u8]) -> Result<u64> {
    left.checked_add(right).ok_or_else(|| {
        RsomicsError::InvalidInput(format!(
            "abundance sum overflows u64 for record {:?}",
            String::from_utf8_lossy(id)
        ))
    })
}

fn validate_range<T>(name: &str, minimum: T, maximum: T) -> Result<()>
where
    T: Ord + std::fmt::Display,
{
    if minimum > maximum {
        return Err(RsomicsError::ConfigError(format!(
            "minimum {name} {minimum} exceeds maximum {name} {maximum}"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn direct_library_records_cannot_bypass_fasta_validation() {
        let record = OwnedRecord {
            id: b"x;size=1".to_vec(),
            seq: Vec::new(),
            qual: None,
        };
        assert!(prepare_record(record, 0, CompatibilityProfile::Strict).is_err());
    }
}
