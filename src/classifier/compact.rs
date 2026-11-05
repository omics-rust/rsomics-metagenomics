use rsomics_common::Result;

use super::invalid;
use super::minimizer::hash;
use super::taxonomy::Taxonomy;

pub(super) struct CompactTable {
    pub cells: Vec<u32>,
    pub value_bits: u8,
    pub occupied: usize,
}

impl CompactTable {
    pub fn new(capacity: usize, bits: u8) -> Result<Self> {
        if capacity == 0 || !(1..=31).contains(&bits) {
            return Err(invalid("invalid compact table capacity or value width"));
        }
        Ok(Self {
            cells: vec![0; capacity],
            value_bits: bits,
            occupied: 0,
        })
    }

    pub fn value(&self, cell: u32) -> usize {
        (cell & ((1u32 << self.value_bits) - 1)) as usize
    }

    pub fn get(&self, key: u64) -> usize {
        let hash = hash(key);
        let fingerprint = (hash >> (32 + self.value_bits)) as u32;
        let mut index = (hash % self.cells.len() as u64) as usize;
        for _ in 0..self.cells.len() {
            let cell = self.cells[index];
            if self.value(cell) == 0 {
                return 0;
            }
            if cell >> self.value_bits == fingerprint {
                return self.value(cell);
            }
            index += 1;
            if index == self.cells.len() {
                index = 0;
            }
        }
        0
    }

    pub fn insert(&mut self, key: u64, taxon: usize, taxonomy: &Taxonomy) -> Result<()> {
        let hash = hash(key);
        let fingerprint = (hash >> (32 + self.value_bits)) as u32;
        let mut index = (hash % self.cells.len() as u64) as usize;
        for _ in 0..self.cells.len() {
            let cell = self.cells[index];
            let value = self.value(cell);
            if value == 0 || cell >> self.value_bits == fingerprint {
                let value = taxonomy.lca(value, taxon);
                self.cells[index] = (fingerprint << self.value_bits) | value as u32;
                self.occupied += usize::from(cell == 0);
                return Ok(());
            }
            index += 1;
            if index == self.cells.len() {
                index = 0;
            }
        }
        Err(invalid("compact hash table capacity exceeded"))
    }
}
