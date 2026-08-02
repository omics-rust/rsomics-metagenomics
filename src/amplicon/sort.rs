use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::OwnedRecord;

use super::{Amplicon, CompatibilityProfile, prepare_record, require_abundance, validate_range};

#[derive(Clone, Debug)]
pub struct SortOptions {
    pub profile: CompatibilityProfile,
    pub minimum_length: usize,
    pub maximum_length: usize,
    pub minimum_abundance: u64,
    pub maximum_abundance: u64,
    pub top: Option<usize>,
}

impl Default for SortOptions {
    fn default() -> Self {
        Self {
            profile: CompatibilityProfile::Strict,
            minimum_length: 1,
            maximum_length: 50_000,
            minimum_abundance: 1,
            maximum_abundance: u64::MAX,
            top: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct AbundanceSort {
    records: Vec<Amplicon>,
    input_records: usize,
    discarded_length: usize,
    discarded_abundance: usize,
}

impl AbundanceSort {
    #[must_use]
    pub fn records(&self) -> &[Amplicon] {
        &self.records
    }

    #[must_use]
    pub fn input_records(&self) -> usize {
        self.input_records
    }

    #[must_use]
    pub fn discarded_length(&self) -> usize {
        self.discarded_length
    }

    #[must_use]
    pub fn discarded_abundance(&self) -> usize {
        self.discarded_abundance
    }
}

pub fn sort_by_abundance<I>(records: I, options: &SortOptions) -> Result<AbundanceSort>
where
    I: IntoIterator<Item = OwnedRecord>,
{
    validate_range("length", options.minimum_length, options.maximum_length)?;
    validate_range(
        "abundance",
        options.minimum_abundance,
        options.maximum_abundance,
    )?;
    if matches!(options.top, Some(0)) {
        return Err(RsomicsError::ConfigError(
            "top must be greater than zero".into(),
        ));
    }

    let mut output = Vec::new();
    let mut input_records = 0usize;
    let mut discarded_length = 0usize;
    let mut discarded_abundance = 0usize;
    for (input_index, record) in records.into_iter().enumerate() {
        input_records += 1;
        let prepared = prepare_record(record, input_index, options.profile)?;
        let abundance = match options.profile {
            CompatibilityProfile::Strict => require_abundance(&prepared)?,
            CompatibilityProfile::Vsearch => prepared.observed_abundance.unwrap_or(1),
        };
        if prepared.sequence.len() < options.minimum_length
            || prepared.sequence.len() > options.maximum_length
        {
            discarded_length += 1;
            continue;
        }
        if abundance < options.minimum_abundance || abundance > options.maximum_abundance {
            discarded_abundance += 1;
            continue;
        }
        let sort_key = if options.profile == CompatibilityProfile::Vsearch {
            prepared.raw_id
        } else {
            prepared.stripped_id.clone()
        };
        output.push(Amplicon {
            id: prepared.stripped_id,
            sequence: prepared.sequence,
            abundance,
            sort_key,
            input_index,
            hash_next: None,
        });
    }
    output.sort_by(|left, right| {
        right
            .abundance
            .cmp(&left.abundance)
            .then_with(|| left.sort_key.cmp(&right.sort_key))
            .then_with(|| left.input_index.cmp(&right.input_index))
    });
    if let Some(top) = options.top {
        output.truncate(top);
    }
    Ok(AbundanceSort {
        records: output,
        input_records,
        discarded_length,
        discarded_abundance,
    })
}
