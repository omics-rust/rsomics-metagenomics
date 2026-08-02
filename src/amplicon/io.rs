use std::io::{BufWriter, Cursor, Read, Write};

use rsomics_common::{Context, Result, RsomicsError};
use rsomics_seqio::{Format, OwnedRecord};

pub fn visit_fasta<R, F>(mut source: R, allow_empty: bool, mut visitor: F) -> Result<usize>
where
    R: Read,
    F: FnMut(OwnedRecord) -> Result<()>,
{
    let mut first = [0u8; 1];
    let count = source
        .read(&mut first)
        .rs_context("reading sequence input")?;
    if count == 0 {
        return if allow_empty {
            Ok(0)
        } else {
            Err(RsomicsError::InvalidInput("FASTA input is empty".into()))
        };
    }
    let replay = Cursor::new(first).chain(source);
    let mut reader = rsomics_seqio::open_reader(replay)?;
    if reader.format() != Format::Fasta {
        return Err(RsomicsError::InvalidInput(
            "amplicon operations require FASTA input".into(),
        ));
    }
    let mut records = 0usize;
    while let Some(record) = reader.read_record()? {
        visitor(record.to_owned())?;
        records = records.checked_add(1).ok_or_else(|| {
            RsomicsError::InvalidInput("FASTA record count overflows usize".into())
        })?;
    }
    Ok(records)
}

pub fn read_fasta<R: Read>(source: R, allow_empty: bool) -> Result<Vec<OwnedRecord>> {
    let mut records = Vec::new();
    visit_fasta(source, allow_empty, |record| {
        records.push(record);
        Ok(())
    })?;
    Ok(records)
}

pub struct FastaWriter<W: Write> {
    inner: BufWriter<W>,
    line_width: usize,
}

impl<W: Write> FastaWriter<W> {
    #[must_use]
    pub fn new(inner: W, line_width: usize) -> Self {
        Self {
            inner: BufWriter::with_capacity(256 * 1024, inner),
            line_width,
        }
    }

    pub fn write_record(&mut self, id: &[u8], sequence: &[u8]) -> Result<()> {
        self.inner.write_all(b">")?;
        self.inner.write_all(id)?;
        self.inner.write_all(b"\n")?;
        if self.line_width == 0 {
            self.inner.write_all(sequence)?;
            self.inner.write_all(b"\n")?;
        } else {
            for line in sequence.chunks(self.line_width) {
                self.inner.write_all(line)?;
                self.inner.write_all(b"\n")?;
            }
        }
        Ok(())
    }

    pub fn finish(mut self) -> Result<W> {
        self.inner.flush()?;
        self.inner
            .into_inner()
            .map_err(|error| RsomicsError::Io(error.into_error()))
    }
}

pub fn fasta_record_size(id_len: usize, sequence_len: usize, line_width: usize) -> Result<u64> {
    let sequence_lines = if line_width == 0 {
        1usize
    } else {
        sequence_len.div_ceil(line_width)
    };
    let bytes = 1usize
        .checked_add(id_len)
        .and_then(|value| value.checked_add(1))
        .and_then(|value| value.checked_add(sequence_len))
        .and_then(|value| value.checked_add(sequence_lines))
        .ok_or_else(|| RsomicsError::InvalidInput("FASTA record size overflows usize".into()))?;
    u64::try_from(bytes)
        .map_err(|_| RsomicsError::InvalidInput("FASTA record size overflows u64".into()))
}
