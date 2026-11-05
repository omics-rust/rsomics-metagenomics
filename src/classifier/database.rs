use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use sha2::{Digest, Sha256};

use super::compact::CompactTable;
use super::taxonomy::Taxonomy;
use super::{DatabaseManifest, invalid};

const HEADER_SIZE: usize = 128;
const MAGIC: &[u8; 8] = b"RSOMDB01";

pub(super) struct Database {
    pub manifest: DatabaseManifest,
    pub taxonomy: Taxonomy,
    pub table: CompactTable,
}

pub(super) fn save(
    path: &Path,
    manifest: &DatabaseManifest,
    taxonomy: &Taxonomy,
    table: &CompactTable,
) -> Result<()> {
    let meta = serde_json::to_vec(manifest).map_err(|error| invalid(error.to_string()))?;
    let nodes = serde_json::to_vec(taxonomy).map_err(|error| invalid(error.to_string()))?;
    let mut staged = super::staged_output(path)?;
    let mut header = [0u8; HEADER_SIZE];
    header[..8].copy_from_slice(MAGIC);
    header[8..16].copy_from_slice(&(meta.len() as u64).to_le_bytes());
    header[16..24].copy_from_slice(&(nodes.len() as u64).to_le_bytes());
    header[24..32].copy_from_slice(&(table.cells.len() as u64).to_le_bytes());
    header[32..64].copy_from_slice(&Sha256::digest(&meta));
    header[64..96].copy_from_slice(&Sha256::digest(&nodes));
    let mut digest = Sha256::new();
    {
        let mut writer = BufWriter::with_capacity(1024 * 1024, staged.as_file_mut());
        writer.write_all(&header)?;
        writer.write_all(&meta)?;
        writer.write_all(&nodes)?;
        let mut bytes = [0u8; 64 * 1024];
        for cells in table.cells.chunks(bytes.len() / 4) {
            for (slot, &cell) in bytes.chunks_exact_mut(4).zip(cells) {
                slot.copy_from_slice(&cell.to_le_bytes());
            }
            let bytes = &bytes[..cells.len() * 4];
            digest.update(bytes);
            writer.write_all(bytes)?;
        }
        writer.flush()?;
    }
    header[96..128].copy_from_slice(&digest.finalize());
    staged.as_file_mut().seek(SeekFrom::Start(0))?;
    staged.as_file_mut().write_all(&header)?;
    staged.as_file_mut().flush()?;
    staged
        .persist_noclobber(path)
        .map_err(|error| RsomicsError::Io(error.error))?;
    Ok(())
}

impl Database {
    pub fn load(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let length = file.metadata()?.len();
        let mut reader = BufReader::new(file);
        let mut header = [0u8; HEADER_SIZE];
        reader.read_exact(&mut header)?;
        if &header[..8] != MAGIC {
            return Err(invalid("unsupported native database format"));
        }
        let number = |offset| {
            u64::from_le_bytes(header[offset..offset + 8].try_into().expect("header field"))
        };
        let (meta_len, node_len, capacity) = (number(8), number(16), number(24));
        let expected = capacity
            .checked_mul(4)
            .and_then(|n| n.checked_add(meta_len))
            .and_then(|n| n.checked_add(node_len))
            .and_then(|n| n.checked_add(HEADER_SIZE as u64));
        if expected != Some(length) {
            return Err(invalid("native database section lengths differ"));
        }
        let mut meta =
            vec![0; usize::try_from(meta_len).map_err(|_| invalid("metadata too large"))?];
        reader.read_exact(&mut meta)?;
        let mut nodes =
            vec![0; usize::try_from(node_len).map_err(|_| invalid("taxonomy too large"))?];
        reader.read_exact(&mut nodes)?;
        if Sha256::digest(&meta)[..] != header[32..64]
            || Sha256::digest(&nodes)[..] != header[64..96]
        {
            return Err(invalid(
                "native database metadata/taxonomy integrity check failed",
            ));
        }
        let manifest: DatabaseManifest =
            serde_json::from_slice(&meta).map_err(|error| invalid(error.to_string()))?;
        if manifest.format_version != 1
            || manifest.profile != "kraken2-2.17.1-linear-v1"
            || manifest.capacity as u64 != capacity
            || manifest.block_size == 0
            || manifest.subblock_size == 0
            || manifest.subblock_size > manifest.block_size
        {
            return Err(invalid("incompatible database parameters"));
        }
        manifest.minimizers.validate()?;
        let taxonomy: Taxonomy =
            serde_json::from_slice(&nodes).map_err(|error| invalid(error.to_string()))?;
        let taxonomy = Taxonomy::from_nodes(taxonomy.nodes)?;
        let mut table = CompactTable::new(
            usize::try_from(capacity).map_err(|_| invalid("capacity too large"))?,
            manifest.value_bits,
        )?;
        if taxonomy.nodes.len() > (1usize << manifest.value_bits) {
            return Err(invalid("taxonomy exceeds compact value width"));
        }
        let mut digest = Sha256::new();
        let mut bytes = [0u8; 64 * 1024];
        for cells in table.cells.chunks_mut(bytes.len() / 4) {
            let bytes = &mut bytes[..cells.len() * 4];
            reader.read_exact(bytes)?;
            digest.update(&*bytes);
            for (cell, raw) in cells.iter_mut().zip(bytes.chunks_exact(4)) {
                *cell = u32::from_le_bytes(raw.try_into().expect("four-byte cell"));
            }
        }
        if digest.finalize()[..] != header[96..128] {
            return Err(invalid("native database table integrity check failed"));
        }
        for &cell in &table.cells {
            let value = table.value(cell);
            if value >= taxonomy.nodes.len() || (value == 0 && cell != 0) {
                return Err(invalid("invalid compact taxonomy value"));
            }
            table.occupied += usize::from(value != 0);
        }
        if table.occupied != manifest.occupied {
            return Err(invalid("database occupancy differs from manifest"));
        }
        Ok(Self {
            manifest,
            taxonomy,
            table,
        })
    }
}
