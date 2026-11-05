use std::collections::VecDeque;

use rsomics_common::Result;
use serde::{Deserialize, Serialize};

use super::invalid;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MinimizerOptions {
    pub k: usize,
    pub length: usize,
    pub spaces: usize,
}

impl Default for MinimizerOptions {
    fn default() -> Self {
        Self {
            k: 35,
            length: 31,
            spaces: 7,
        }
    }
}

impl MinimizerOptions {
    pub(super) fn validate(&self) -> Result<()> {
        if self.length == 0
            || self.length > 31
            || self.k < self.length
            || self.spaces > self.length / 4
        {
            return Err(invalid(
                "require 1 <= minimizer length <= min(k,31) and spaces <= length/4",
            ));
        }
        Ok(())
    }
}

pub(super) fn hash(mut value: u64) -> u64 {
    value ^= value >> 33;
    value = value.wrapping_mul(0xff51afd7ed558ccd);
    value ^= value >> 33;
    value = value.wrapping_mul(0xc4ceb9fe1a85ec53);
    value ^ (value >> 33)
}

pub(super) struct Hit {
    pub key: u64,
    pub ambiguous: bool,
}

pub(super) struct Scanner<'a> {
    sequence: &'a [u8],
    k: usize,
    length: usize,
    start: usize,
    finish: usize,
    position: usize,
    mask: u64,
    spaced: u64,
    toggle: u64,
    lmer: u64,
    ambiguity: u64,
    loaded: usize,
    queue_position: usize,
    queue: VecDeque<(u64, usize)>,
    last: u64,
}

impl<'a> Scanner<'a> {
    pub fn new(
        sequence: &'a [u8],
        options: &MinimizerOptions,
        start: usize,
        finish: usize,
    ) -> Self {
        let mask = (1u64 << (2 * options.length)) - 1;
        let mut spaced = mask;
        for index in 0..options.spaces {
            spaced &= !(3u64 << (2 * (index * 2 + 1)));
        }
        Self {
            sequence,
            k: options.k,
            length: options.length,
            start,
            finish: finish.min(sequence.len()),
            position: start,
            mask,
            spaced,
            toggle: 0xe37e28c4271b5a2d & mask,
            lmer: 0,
            ambiguity: 0,
            loaded: 0,
            queue_position: 0,
            queue: VecDeque::with_capacity((options.k - options.length + 1).min(sequence.len())),
            last: u64::MAX,
        }
    }

    fn hit(&self) -> Hit {
        Hit {
            key: self.last,
            ambiguous: self.queue_position < self.k - self.length || self.ambiguity != 0,
        }
    }
}

impl Iterator for Scanner<'_> {
    type Item = Hit;

    fn next(&mut self) -> Option<Hit> {
        if self.position >= self.finish {
            return None;
        }
        loop {
            if self.loaded == self.length {
                self.loaded -= 1;
            }
            while self.loaded < self.length && self.position < self.finish {
                self.loaded += 1;
                self.lmer <<= 2;
                self.ambiguity <<= 2;
                let base = self.sequence[self.position];
                self.position += 1;
                let code = match base {
                    b'A' | b'a' => Some(0),
                    b'C' | b'c' => Some(1),
                    b'G' | b'g' => Some(2),
                    b'T' | b't' => Some(3),
                    _ => None,
                };
                if let Some(code) = code {
                    self.lmer |= code;
                } else {
                    self.queue.clear();
                    self.queue_position = 0;
                    self.lmer = 0;
                    self.loaded = 0;
                    self.ambiguity |= 3;
                }
                self.lmer &= self.mask;
                self.ambiguity &= self.mask;
                if self.position - self.start >= self.k && self.loaded < self.length {
                    return Some(self.hit());
                }
            }
            if self.loaded < self.length {
                return None;
            }
            let reverse = self.lmer.reverse_bits();
            let reverse =
                ((reverse & 0xaaaa_aaaa_aaaa_aaaa) >> 1) | ((reverse & 0x5555_5555_5555_5555) << 1);
            let reverse = ((!reverse) >> (64 - 2 * self.length)) & self.mask;
            let key = self.lmer.min(reverse) & self.spaced;
            let candidate = key ^ self.toggle;
            if self.k == self.length {
                self.last = key;
                return Some(self.hit());
            }
            while self
                .queue
                .back()
                .is_some_and(|&(value, _)| value > candidate)
            {
                self.queue.pop_back();
            }
            self.queue.push_back((candidate, self.queue_position));
            let difference = self.k - self.length;
            if self.queue_position >= difference
                && self
                    .queue
                    .front()
                    .is_some_and(|&(_, position)| position < self.queue_position - difference)
            {
                self.queue.pop_front();
            }
            self.queue_position += 1;
            // Kraken's interval scanner uses the absolute sequence coordinate here.
            if self.position >= self.k {
                break;
            }
        }
        self.last = self.queue.front().expect("complete minimizer window").0 ^ self.toggle;
        Some(self.hit())
    }
}
