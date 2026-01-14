use crate::data::{ByteIndexMap, matrices::WeightMatrix};

/// Parses a matrix from a file.
///
/// Comment lines beginning with `#` are ignored at the beginning of the file.
/// The alphabet is parsed from the first non-comment line, and must have at
/// least the residues specified in `byte_index_map`. Row labels are optionally
/// allowed.
#[must_use]
pub const fn parse_matrix<const S: usize>(
    bytes: &[u8], byte_index_map: &'static ByteIndexMap<S>,
) -> WeightMatrix<'static, i8, S> {
    MatParser::new(bytes).parse_bla_matrix(byte_index_map)
}

/// A struct for parsing a protein substitution matrix at compile time.
///
/// This tracks the original data and the current position within that data.
/// With the exception of [`get`] and [`advance`], any methods that return an
/// `Option` or `bool` to indicate the potential of failure may invalidate the
/// position stored in the [`MatParser`].
///
/// [`get`]: MatParser::get
/// [`advance`]: MatParser::advance
pub(crate) struct MatParser<'a> {
    bytes:    &'a [u8],
    position: usize,
}

impl<'a> MatParser<'a> {
    /// Creates a new [`MatParser`] from the provided bytes.
    pub(crate) const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, position: 0 }
    }

    /// Returns the length of the bytes in the underlying data.
    const fn len(&self) -> usize {
        self.bytes.len()
    }

    /// Retrieves the byte at the current position if it exists, without
    /// advancing to the next position.
    const fn get(&self) -> Option<u8> {
        if self.position < self.len() {
            Some(self.bytes[self.position])
        } else {
            None
        }
    }

    /// Retrieves the byte at the current position if it exists, and if it does,
    /// advance to the next position.
    const fn advance(&mut self) -> Option<u8> {
        let Some(out) = self.get() else {
            return None;
        };
        self.position += 1;
        Some(out)
    }

    /// Advances the position past any spaces or tabs.
    const fn skip_spaces(&mut self) {
        while let Some(byte) = self.get()
            && (byte == b' ' || byte == b'\t')
        {
            self.position += 1;
        }
    }

    /// Advances the position past the next line break, either `\n` or `\r\n`.
    /// If neither of these are found, `false` is returned.
    const fn parse_linebreak(&mut self) -> bool {
        let Some(byte) = self.advance() else {
            return false;
        };

        if byte == b'\n' {
            return true;
        }

        if byte == b'\r' {
            let Some(byte) = self.advance() else {
                return false;
            };
            if byte == b'\n' {
                return true;
            }
        }

        false
    }

    /// Advances the position past the next line break or to the end of the
    /// data. If an invalid line break is encountered (`\r` without a subsequent
    /// `\n`), `false` is returned.
    const fn parse_rest_of_line(&mut self) -> bool {
        while let Some(byte) = self.get() {
            if byte == b'\r' || byte == b'\n' {
                return self.parse_linebreak();
            }
            self.position += 1;
        }
        true
    }

    /// Parses the next byte in the data as a digit. `None` is returned if this
    /// fails.
    const fn parse_digit(&mut self) -> Option<u8> {
        let Some(byte) = self.advance() else {
            return None;
        };
        if !byte.is_ascii_digit() {
            return None;
        }
        Some(byte - b'0')
    }

    /// Parses the next bytes as an `i8`.
    #[allow(clippy::cast_possible_wrap)]
    const fn parse_i8(&mut self) -> Option<i8> {
        let Some(byte) = self.get() else {
            return None;
        };
        let is_negative = byte == b'-';

        if is_negative {
            self.advance();
        }

        let Some(mut magnitude) = self.parse_digit() else {
            return None;
        };

        while let Some(byte) = self.get() {
            if byte.is_ascii_whitespace() {
                break;
            }
            let Some(new_digit) = self.parse_digit() else {
                return None;
            };

            let Some(magnitude_t_10) = magnitude.checked_mul(10u8) else {
                return None;
            };
            let Some(new_magnitude) = magnitude_t_10.checked_add(new_digit) else {
                return None;
            };
            magnitude = new_magnitude;
        }

        if is_negative {
            0i8.checked_sub_unsigned(magnitude)
        } else {
            0i8.checked_add_unsigned(magnitude)
        }
    }

    /// Parses a row of whitespace-separated `i8` scores, returning `None` if
    /// insufficiently many values are found or a parsing error occurs.
    ///
    /// ## Panics
    ///
    /// If `row_label` is provided and the first non-whitespace character is not
    /// `row_label`, then this will panic.
    const fn parse_row<const S: usize>(&mut self, row_label: Option<u8>) -> Option<[i8; S]> {
        self.skip_spaces();

        if let Some(row_label) = row_label {
            let Some(label) = self.advance() else {
                return None;
            };
            assert!(
                label == row_label,
                "Invalid row label found! If row labels are provided, ensure the row labels and column labels are in the same order"
            );
            self.skip_spaces();
        }

        let mut out = [0; S];
        let mut i = 0;

        while i < S {
            let Some(new_element) = self.parse_i8() else {
                panic!("Failed to parse a number in a row of the weight matrix")
            };
            self.skip_spaces();

            out[i] = new_element;
            i += 1;
        }

        assert!(i >= S, "Too few values found in a row of the weight matrix");

        Some(out)
    }

    /// Parses a matrix of `i8` scores where the columns are
    /// whitespace-separated and the rows each are on a new line. `None` is
    /// returned if insufficiently many rows or columns are found, or if a
    /// parsing error occurs.
    ///
    /// If `row_labels` is provided, then it is expected that the matrix will
    /// contain row labels which are in the same order as the column labels.
    const fn parse_mat<const S: usize>(&mut self, row_labels: Option<[u8; S]>) -> Option<[[i8; S]; S]> {
        assert!(S > 0, "Cannot have an empty alphabet");

        let mut out = [[0; S]; S];

        let row_label = match row_labels {
            Some(row_labels) => Some(row_labels[0]),
            None => None,
        };
        let Some(row) = self.parse_row(row_label) else {
            panic!("Could not parse first row in weight matrix");
        };
        out[0] = row;

        let mut i = 1;
        while i < S && self.position < self.len() {
            if !self.parse_linebreak() {
                return None;
            }

            let row_label = match row_labels {
                Some(row_labels) => Some(row_labels[i]),
                None => None,
            };

            let Some(row) = self.parse_row(row_label) else {
                return None;
            };
            out[i] = row;
            i += 1;
        }

        if i < S {
            return None;
        }

        Some(out)
    }

    /// Parses the alphabet line of a `.bla` file (whitespace-separated
    /// characters).
    ///
    /// ## Panics
    ///
    /// If too many or too few characters are found, or the characters are not
    /// graphic ASCII, this panics.
    const fn parse_alphabet_line<const S: usize>(&mut self) -> [u8; S] {
        self.skip_spaces();

        let mut out = [0; S];
        let mut i = 0;

        while i < S {
            let Some(new_element) = self.advance() else {
                panic!("Failed to parse the alphabet line: too few characters found")
            };
            assert!(
                !(new_element == b'\r' || new_element == b'\n'),
                "Failed to parse the alphabet line: too few characters found"
            );
            assert!(new_element.is_ascii_graphic(), "Invalid character found!");
            self.skip_spaces();

            out[i] = new_element;
            i += 1;
        }

        assert!(
            self.parse_linebreak(),
            "Failed to parse the line break after the alphabet line"
        );

        out
    }

    /// Before parsing a matrix, check whether it contains row labels.
    const fn has_row_labels(&mut self) -> bool {
        let original_position = self.position;
        self.skip_spaces();
        let Some(mut byte) = self.advance() else {
            self.position = original_position;
            return false;
        };
        if byte == b'-' {
            let Some(new_byte) = self.advance() else {
                self.position = original_position;
                return true;
            };
            byte = new_byte;
        }
        self.position = original_position;
        !byte.is_ascii_digit()
    }

    /// Parses a matrix from a `.bla` file. Comment lines are ignored. The
    /// alphabet is parsed from the first non-comment line. Row labels are
    /// optionally allowed.
    pub(crate) const fn parse_bla_matrix<const S: usize>(
        &mut self, byte_index_map: &'static ByteIndexMap<S>,
    ) -> WeightMatrix<'static, i8, S> {
        assert!(S > 0, "Cannot have an empty alphabet");

        self.skip_spaces();
        while let Some(byte) = self.get() {
            if byte == b'#' {
                assert!(self.parse_rest_of_line(), "Invalid line break in file");
            } else {
                break;
            }
        }
        let alphabet = self.parse_alphabet_line::<S>();

        // Case sensitivity and catch_all do not matter, since we will be
        // subsetting this
        let file_byte_index_map = ByteIndexMap::new(alphabet, alphabet[0]);

        let has_row_labels = self.has_row_labels();
        let row_labels = if has_row_labels { Some(alphabet) } else { None };

        let Some(weights) = self.parse_mat::<S>(row_labels) else {
            panic!("Failed to parse weight matrix");
        };

        let file_weight_matrix = WeightMatrix::new_custom(&file_byte_index_map, weights);
        file_weight_matrix.get_subset(byte_index_map)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse_i8() {
        let mut m = MatParser::new(b"-128");
        assert_eq!(Some(-128i8), m.parse_i8());
    }
}
