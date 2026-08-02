use rsomics_common::{Result, RsomicsError};

use super::CompatibilityProfile;

pub(super) struct ParsedHeader {
    pub raw: Vec<u8>,
    pub abundance: Option<u64>,
    selected: Option<(usize, usize)>,
}

impl ParsedHeader {
    pub fn stripped(&self) -> Vec<u8> {
        self.selected.map_or_else(
            || self.raw.clone(),
            |range| remove_attribute(&self.raw, range),
        )
    }
}

pub(super) fn parse_header(mut id: Vec<u8>, profile: CompatibilityProfile) -> Result<ParsedHeader> {
    if profile == CompatibilityProfile::Vsearch
        && let Some(end) = id.iter().position(|byte| byte.is_ascii_whitespace())
    {
        id.truncate(end);
    }
    let raw = id;

    let mut abundance = None;
    let mut selected = None;
    let mut segment_start = 0usize;
    loop {
        let segment_end = raw[segment_start..]
            .iter()
            .position(|byte| *byte == b';')
            .map_or(raw.len(), |offset| segment_start + offset);
        let segment = &raw[segment_start..segment_end];
        if let Some(digits) = segment.strip_prefix(b"size=") {
            match profile {
                CompatibilityProfile::Strict => {
                    let value = parse_decimal(digits, u64::MAX, &raw)?;
                    if abundance.replace(value).is_some() {
                        return Err(invalid_header(&raw, "duplicate size abundance"));
                    }
                    selected = Some((segment_start, segment_end));
                }
                CompatibilityProfile::Vsearch => {
                    if abundance.is_none()
                        && !digits.is_empty()
                        && digits.iter().all(u8::is_ascii_digit)
                    {
                        let value = parse_decimal(digits, i64::MAX as u64, &raw)?;
                        abundance = Some(value);
                        selected = Some((segment_start, segment_end));
                    }
                }
            }
        }
        if segment_end == raw.len() {
            break;
        }
        segment_start = segment_end + 1;
    }

    Ok(ParsedHeader {
        raw,
        abundance,
        selected,
    })
}

fn parse_decimal(digits: &[u8], limit: u64, header: &[u8]) -> Result<u64> {
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return Err(invalid_header(header, "malformed size abundance"));
    }
    let mut value = 0u64;
    for digit in digits {
        value = value
            .checked_mul(10)
            .and_then(|value| value.checked_add(u64::from(digit - b'0')))
            .filter(|value| *value <= limit)
            .ok_or_else(|| invalid_header(header, "size abundance is out of range"))?;
    }
    if value == 0 {
        return Err(invalid_header(
            header,
            "size abundance must be greater than zero",
        ));
    }
    Ok(value)
}

fn remove_attribute(header: &[u8], (start, end): (usize, usize)) -> Vec<u8> {
    let (remove_start, remove_end) = if start > 0 {
        (start - 1, end)
    } else if end < header.len() {
        (0, end + 1)
    } else {
        (0, end)
    };
    let mut stripped = Vec::with_capacity(header.len() - (remove_end - remove_start));
    stripped.extend_from_slice(&header[..remove_start]);
    stripped.extend_from_slice(&header[remove_end..]);
    stripped
}

fn invalid_header(header: &[u8], reason: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "{reason} in record {:?}",
        String::from_utf8_lossy(header)
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strict_rejects_duplicate_and_malformed_values() {
        assert!(parse_header(b"x;size=2;size=3".to_vec(), CompatibilityProfile::Strict).is_err());
        assert!(parse_header(b"x;size=no".to_vec(), CompatibilityProfile::Strict).is_err());
        assert!(parse_header(b"x;size=0".to_vec(), CompatibilityProfile::Strict).is_err());
        assert!(
            parse_header(
                b"x;size=18446744073709551616".to_vec(),
                CompatibilityProfile::Strict,
            )
            .is_err()
        );
        assert!(
            parse_header(
                b"x;size=9223372036854775808".to_vec(),
                CompatibilityProfile::Vsearch,
            )
            .is_err()
        );
    }

    #[test]
    fn vsearch_selects_and_removes_only_the_first_valid_value() {
        let parsed = parse_header(
            b"x;size=no;size=2;size=3".to_vec(),
            CompatibilityProfile::Vsearch,
        )
        .unwrap();
        assert_eq!(parsed.abundance, Some(2));
        assert_eq!(parsed.stripped(), b"x;size=no;size=3");
    }

    #[test]
    fn leading_and_middle_attributes_strip_one_separator() {
        let leading =
            parse_header(b"size=2;sample=A".to_vec(), CompatibilityProfile::Strict).unwrap();
        let middle =
            parse_header(b"x;size=2;sample=A".to_vec(), CompatibilityProfile::Strict).unwrap();
        assert_eq!(leading.stripped(), b"sample=A");
        assert_eq!(middle.stripped(), b"x;sample=A");
    }
}
