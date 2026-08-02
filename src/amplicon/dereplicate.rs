use ahash::AHashMap;
use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::OwnedRecord;
use std::io::Read;

use super::{
    Amplicon, CompatibilityProfile, PreparedRecord, checked_sum, prepare_record, require_abundance,
    validate_range, visit_fasta,
};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum DereplicateMode {
    #[default]
    FullLength,
    Prefix,
}

#[derive(Clone, Debug)]
pub struct DereplicateOptions {
    pub mode: DereplicateMode,
    pub profile: CompatibilityProfile,
    pub input_abundance: bool,
    pub minimum_length: usize,
    pub maximum_length: usize,
    pub minimum_abundance: u64,
    pub maximum_abundance: u64,
    pub top: Option<usize>,
}

impl Default for DereplicateOptions {
    fn default() -> Self {
        Self {
            mode: DereplicateMode::FullLength,
            profile: CompatibilityProfile::Strict,
            input_abundance: false,
            minimum_length: 32,
            maximum_length: 50_000,
            minimum_abundance: 1,
            maximum_abundance: u64::MAX,
            top: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Dereplication {
    records: Vec<Amplicon>,
    input_records: usize,
    discarded_length: usize,
    clusters_before_filter: usize,
}

impl Dereplication {
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
    pub fn clusters_before_filter(&self) -> usize {
        self.clusters_before_filter
    }
}

pub fn dereplicate<I>(records: I, options: &DereplicateOptions) -> Result<Dereplication>
where
    I: IntoIterator<Item = OwnedRecord>,
{
    validate_options(options)?;
    match options.mode {
        DereplicateMode::FullLength => full_length(records, options),
        DereplicateMode::Prefix => prefix(records, options),
    }
}

pub fn dereplicate_fasta<R: Read>(
    source: R,
    options: &DereplicateOptions,
) -> Result<Dereplication> {
    validate_options(options)?;
    if options.mode == DereplicateMode::Prefix {
        let records = super::read_fasta(source, options.profile == CompatibilityProfile::Vsearch)?;
        return prefix(records, options);
    }
    let mut state = FullLengthState::default();
    visit_fasta(
        source,
        options.profile == CompatibilityProfile::Vsearch,
        |record| state.push(record, options),
    )?;
    state.finish(options)
}

#[derive(Default)]
struct FullLengthState {
    index: AHashMap<u64, usize>,
    clusters: Vec<Amplicon>,
    input_records: usize,
    discarded_length: usize,
}

struct PrefixRecord {
    prepared: PreparedRecord,
    normalized: Option<Vec<u8>>,
    abundance: u64,
    ordering_abundance: u64,
}

impl PrefixRecord {
    fn normalized(&self) -> &[u8] {
        self.normalized
            .as_deref()
            .unwrap_or(&self.prepared.sequence)
    }
}

impl FullLengthState {
    fn push(&mut self, record: OwnedRecord, options: &DereplicateOptions) -> Result<()> {
        let input_index = self.input_records;
        self.input_records = self.input_records.checked_add(1).ok_or_else(|| {
            RsomicsError::InvalidInput("FASTA record count overflows usize".into())
        })?;
        if !length_in_range(record.seq.len(), options) {
            self.discarded_length += 1;
            return Ok(());
        }
        let prepared = prepare_record(record, input_index, options.profile)?;
        let abundance = input_abundance(&prepared, options)?;
        let hash = normalized_hash(&prepared.sequence, options.profile, &prepared.raw_id)?;
        let mut candidate = self.index.get(&hash).copied();
        let mut matching = None;
        while let Some(cluster_index) = candidate {
            let cluster = &self.clusters[cluster_index];
            if normalized_equal(&prepared.sequence, &cluster.sequence) {
                matching = Some(cluster_index);
                break;
            }
            candidate = cluster.hash_next;
        }
        if let Some(cluster_index) = matching {
            let cluster = &mut self.clusters[cluster_index];
            cluster.abundance = checked_sum(cluster.abundance, abundance, cluster.id())?;
        } else {
            let cluster_index = self.clusters.len();
            let hash_next = self.index.insert(hash, cluster_index);
            let mut cluster = to_amplicon(prepared, abundance, options.profile);
            cluster.hash_next = hash_next;
            self.clusters.push(cluster);
        }
        Ok(())
    }

    fn finish(self, options: &DereplicateOptions) -> Result<Dereplication> {
        finish(
            self.clusters,
            self.input_records,
            self.discarded_length,
            options,
        )
    }
}

fn full_length<I>(records: I, options: &DereplicateOptions) -> Result<Dereplication>
where
    I: IntoIterator<Item = OwnedRecord>,
{
    let mut state = FullLengthState::default();
    for record in records {
        state.push(record, options)?;
    }
    state.finish(options)
}

fn prefix<I>(records: I, options: &DereplicateOptions) -> Result<Dereplication>
where
    I: IntoIterator<Item = OwnedRecord>,
{
    let mut input_records = 0usize;
    let mut discarded_length = 0usize;
    let mut raw = Vec::<PrefixRecord>::new();
    for (input_index, record) in records.into_iter().enumerate() {
        input_records += 1;
        if !length_in_range(record.seq.len(), options) {
            discarded_length += 1;
            continue;
        }
        let prepared = prepare_record(record, input_index, options.profile)?;
        let abundance = input_abundance(&prepared, options)?;
        let ordering_abundance = if options.profile == CompatibilityProfile::Vsearch {
            prepared.observed_abundance.unwrap_or(1)
        } else {
            abundance
        };
        let normalized = normalized_copy(&prepared.sequence, options.profile, &prepared.raw_id)?;
        raw.push(PrefixRecord {
            prepared,
            normalized,
            abundance,
            ordering_abundance,
        });
    }

    raw.sort_by(|left, right| {
        left.normalized()
            .len()
            .cmp(&right.normalized().len())
            .then_with(|| right.ordering_abundance.cmp(&left.ordering_abundance))
            .then_with(|| left.prepared.raw_id.cmp(&right.prepared.raw_id))
            .then_with(|| left.prepared.input_index.cmp(&right.prepared.input_index))
    });
    if raw.is_empty() {
        return finish(Vec::new(), input_records, discarded_length, options);
    }

    #[derive(Clone, Default)]
    struct Bucket {
        hash: u64,
        record: Option<usize>,
        abundance: u64,
        deleted: bool,
    }

    let table_size = table_size(raw.len())?;
    let mask = table_size - 1;
    let mut table = vec![Bucket::default(); table_size];
    let shortest = raw[0].normalized().len();

    for record_index in 0..raw.len() {
        let hashes = prefix_hashes(raw[record_index].normalized());
        let length = raw[record_index].normalized().len();
        let exact_hash = hashes[length];
        let mut exact_slot = slot(exact_hash, mask);
        loop {
            let bucket = &table[exact_slot];
            let Some(existing) = bucket.record else {
                break;
            };
            if !bucket.deleted
                && bucket.hash == exact_hash
                && raw[existing].normalized() == raw[record_index].normalized()
            {
                table[exact_slot].abundance = checked_sum(
                    table[exact_slot].abundance,
                    raw[record_index].abundance,
                    &raw[existing].prepared.stripped_id,
                )?;
                break;
            }
            exact_slot = (exact_slot + 1) & mask;
        }
        if table[exact_slot].record.is_some() {
            continue;
        }

        let mut replacement = None;
        for prefix_length in (shortest..length).rev() {
            let hash = hashes[prefix_length];
            let mut candidate_slot = slot(hash, mask);
            loop {
                let bucket = &table[candidate_slot];
                let Some(existing) = bucket.record else {
                    break;
                };
                if !bucket.deleted
                    && bucket.hash == hash
                    && raw[existing].normalized().len() == prefix_length
                    && raw[existing].normalized()
                        == &raw[record_index].normalized()[..prefix_length]
                {
                    replacement = Some(candidate_slot);
                    break;
                }
                candidate_slot = (candidate_slot + 1) & mask;
            }
            if replacement.is_some() {
                break;
            }
        }

        let inherited = if let Some(replaced) = replacement {
            table[replaced].deleted = true;
            table[replaced].abundance
        } else {
            0
        };
        table[exact_slot] = Bucket {
            hash: exact_hash,
            record: Some(record_index),
            abundance: checked_sum(
                inherited,
                raw[record_index].abundance,
                &raw[record_index].prepared.stripped_id,
            )?,
            deleted: false,
        };
    }

    let clusters = table
        .into_iter()
        .filter(|bucket| bucket.record.is_some() && !bucket.deleted)
        .map(|bucket| {
            let record = bucket.record.expect("filtered populated bucket");
            to_amplicon(
                raw[record].prepared.clone(),
                bucket.abundance,
                options.profile,
            )
        })
        .collect();
    finish(clusters, input_records, discarded_length, options)
}

fn finish(
    mut records: Vec<Amplicon>,
    input_records: usize,
    discarded_length: usize,
    options: &DereplicateOptions,
) -> Result<Dereplication> {
    let clusters_before_filter = records.len();
    records.retain(|record| {
        record.abundance >= options.minimum_abundance
            && record.abundance <= options.maximum_abundance
    });
    records.sort_by(|left, right| {
        right
            .abundance
            .cmp(&left.abundance)
            .then_with(|| left.sort_key.cmp(&right.sort_key))
            .then_with(|| left.input_index.cmp(&right.input_index))
    });
    if let Some(top) = options.top {
        records.truncate(top);
    }
    Ok(Dereplication {
        records,
        input_records,
        discarded_length,
        clusters_before_filter,
    })
}

fn input_abundance(record: &PreparedRecord, options: &DereplicateOptions) -> Result<u64> {
    if options.input_abundance {
        match options.profile {
            CompatibilityProfile::Strict => require_abundance(record),
            CompatibilityProfile::Vsearch => Ok(record.observed_abundance.unwrap_or(1)),
        }
    } else {
        Ok(1)
    }
}

fn to_amplicon(record: PreparedRecord, abundance: u64, profile: CompatibilityProfile) -> Amplicon {
    let sort_key = if profile == CompatibilityProfile::Vsearch {
        record.raw_id
    } else {
        record.stripped_id.clone()
    };
    Amplicon {
        id: record.stripped_id,
        sequence: record.sequence,
        abundance,
        sort_key,
        input_index: record.input_index,
        hash_next: None,
    }
}

fn length_in_range(length: usize, options: &DereplicateOptions) -> bool {
    length >= options.minimum_length && length <= options.maximum_length
}

fn validate_options(options: &DereplicateOptions) -> Result<()> {
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
    Ok(())
}

fn normalized_copy(
    sequence: &[u8],
    profile: CompatibilityProfile,
    id: &[u8],
) -> Result<Option<Vec<u8>>> {
    let mut changed = false;
    for byte in sequence {
        changed |= normalize_byte(*byte, profile, id)? != *byte;
    }
    if !changed {
        return Ok(None);
    }
    sequence
        .iter()
        .copied()
        .map(|byte| normalize_byte(byte, profile, id))
        .collect::<Result<Vec<_>>>()
        .map(Some)
}

fn normalized_hash(sequence: &[u8], profile: CompatibilityProfile, id: &[u8]) -> Result<u64> {
    let mut hash = FNV_OFFSET;
    for byte in sequence {
        hash ^= u64::from(normalize_byte(*byte, profile, id)?);
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    Ok(hash)
}

fn normalized_equal(left: &[u8], right: &[u8]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(left, right)| normalize_known(*left) == normalize_known(*right))
}

fn normalize_known(byte: u8) -> u8 {
    match byte.to_ascii_uppercase() {
        b'U' => b'T',
        value @ (b'A' | b'B' | b'C' | b'D' | b'G' | b'H' | b'K' | b'M' | b'N' | b'R' | b'S'
        | b'T' | b'V' | b'W' | b'Y') => value,
        _ => b'N',
    }
}

fn normalize_byte(byte: u8, profile: CompatibilityProfile, id: &[u8]) -> Result<u8> {
    let normalized = match byte.to_ascii_uppercase() {
        b'U' => b'T',
        value @ (b'A' | b'B' | b'C' | b'D' | b'G' | b'H' | b'K' | b'M' | b'N' | b'R' | b'S'
        | b'T' | b'V' | b'W' | b'Y') => value,
        _ if profile == CompatibilityProfile::Vsearch => b'N',
        _ => {
            return Err(RsomicsError::InvalidInput(format!(
                "record {:?} contains a non-IUPAC nucleotide byte {byte}",
                String::from_utf8_lossy(id)
            )));
        }
    };
    Ok(normalized)
}

const FNV_OFFSET: u64 = 14_695_981_039_346_656_037;
const FNV_PRIME: u64 = 1_099_511_628_211;

fn prefix_hashes(sequence: &[u8]) -> Vec<u64> {
    let mut hashes = Vec::with_capacity(sequence.len() + 1);
    let mut hash = FNV_OFFSET;
    hashes.push(hash);
    for byte in sequence {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(FNV_PRIME);
        hashes.push(hash);
    }
    hashes
}

fn table_size(records: usize) -> Result<usize> {
    let required = records
        .checked_mul(3)
        .and_then(|value| value.checked_add(1))
        .map(|value| value / 2)
        .ok_or_else(|| RsomicsError::InvalidInput("prefix table size overflows usize".into()))?;
    required
        .max(1)
        .checked_next_power_of_two()
        .ok_or_else(|| RsomicsError::InvalidInput("prefix table size overflows usize".into()))
}

fn slot(hash: u64, mask: usize) -> usize {
    (hash as usize) & mask
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(id: &str, sequence: &str) -> OwnedRecord {
        OwnedRecord {
            id: id.as_bytes().to_vec(),
            seq: sequence.as_bytes().to_vec(),
            qual: None,
        }
    }

    #[test]
    fn full_length_normalizes_case_and_u() {
        let options = DereplicateOptions {
            minimum_length: 1,
            input_abundance: true,
            ..DereplicateOptions::default()
        };
        let result = dereplicate(
            [
                record("first;size=2", "acgu"),
                record("second;size=3", "ACGT"),
            ],
            &options,
        )
        .unwrap();
        assert_eq!(result.records().len(), 1);
        assert_eq!(result.records()[0].abundance(), 5);
        assert_eq!(result.records()[0].id(), b"first");
        assert_eq!(result.records()[0].sequence(), b"acgu");
    }

    #[test]
    fn prefix_replaces_shorter_representative() {
        let options = DereplicateOptions {
            mode: DereplicateMode::Prefix,
            minimum_length: 1,
            ..DereplicateOptions::default()
        };
        let result = dereplicate(
            [record("short", "ACGT"), record("long", "ACGTAA")],
            &options,
        )
        .unwrap();
        assert_eq!(result.records().len(), 1);
        assert_eq!(result.records()[0].id(), b"long");
        assert_eq!(result.records()[0].abundance(), 2);
    }

    #[test]
    fn abundance_sum_overflow_is_an_error() {
        let options = DereplicateOptions {
            minimum_length: 1,
            input_abundance: true,
            ..DereplicateOptions::default()
        };
        let result = dereplicate(
            [
                record("first;size=18446744073709551615", "ACGT"),
                record("second;size=1", "ACGT"),
            ],
            &options,
        );
        assert!(result.is_err());
    }
}
