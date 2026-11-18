use crate::search::ByteSubstring;
use std::io::{BufRead, BufReader, ErrorKind, Read};

/// A wrapper around [`BufReader`] that stops yielding bytes permanently after
/// EOF (i.e., it fuses after the inner [`BufReader`] stops yielding bytes).
#[derive(Debug)]
pub struct FusedBufReader<R> {
    /// The inner [`BufReader`], set to `None` after EOF is reached.
    inner: Option<BufReader<R>>,
}

impl<R: Read> FusedBufReader<R> {
    /// Creates a new [`FusedBufReader`] around something implementing [`Read`].
    #[allow(dead_code)]
    pub fn new(readable: R) -> Self {
        Self::from_bufreader(BufReader::new(readable))
    }

    /// Creates a new [`FusedBufReader`] from an existing [`BufReader`].
    #[allow(dead_code)]
    pub fn from_bufreader(bufreader: BufReader<R>) -> Self {
        FusedBufReader { inner: Some(bufreader) }
    }

    /// Returns whether the [`BufReader`] has reached EOF and will therefore
    /// read no more bytes.
    #[allow(dead_code)]
    pub fn is_at_eof(&self) -> bool {
        self.inner.is_none()
    }
}

impl<R: Read> Read for FusedBufReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let Some(inner) = &mut self.inner else {
            return Ok(0);
        };

        let bytes = inner.read(buf)?;
        if bytes == 0 && !buf.is_empty() {
            self.inner = None;
        }
        Ok(bytes)
    }

    fn read_to_end(&mut self, buf: &mut Vec<u8>) -> std::io::Result<usize> {
        let Some(inner) = &mut self.inner else {
            return Ok(0);
        };

        let bytes = inner.read_to_end(buf)?;
        self.inner = None;
        Ok(bytes)
    }

    fn read_to_string(&mut self, buf: &mut String) -> std::io::Result<usize> {
        let Some(inner) = &mut self.inner else {
            return Ok(0);
        };

        let bytes = inner.read_to_string(buf)?;
        self.inner = None;
        Ok(bytes)
    }
}

impl<R: Read> BufRead for FusedBufReader<R> {
    #[inline]
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        let Some(inner) = &mut self.inner else {
            return Ok(&[]);
        };

        let out = inner.fill_buf()?;
        if out.is_empty() {
            self.inner = None;
            Ok(&[])
        } else {
            Ok(out)
        }
    }

    #[inline]
    fn consume(&mut self, amount: usize) {
        if let Some(inner) = &mut self.inner {
            inner.consume(amount);
        }
    }
}

/// An extension trait for [`BufRead`] offering a retrying version of
/// [`BufRead::fill_buf`].
pub trait BufReadExtension: BufRead {
    /// A version of [`BufRead::fill_buf`] that retries upon
    /// [`ErrorKind::Interrupted`].
    ///
    /// Many higher-level [`BufRead`] methods retry on this error like
    /// [`BufRead::read_until`], so this method allows implementing similar
    /// higher-level functionality.
    fn retrying_fill_buf(&mut self) -> std::io::Result<&[u8]> {
        loop {
            match self.fill_buf() {
                Ok(buf) => return Ok(buf),
                Err(err) if err.kind() == ErrorKind::Interrupted => {}
                Err(err) => return Err(err),
            }
        }
    }

    #[inline]
    fn read_line_owned(&mut self, capacity: usize) -> std::io::Result<Vec<u8>> {
        let mut filled = self.retrying_fill_buf()?;

        // If the line is entirely within the buffer, then use to_vec for
        // exact-size allocation
        let mut line = if let Some(i) = filled.find_byte(b'\n') {
            let line = filled[..i].to_vec();
            self.consume(i + 1);
            line
        } else {
            let mut line = Vec::with_capacity(capacity.max(filled.len()));

            loop {
                if filled.is_empty() {
                    return Ok(line);
                }
                line.extend_from_slice(filled);
                let used = filled.len();
                self.consume(used);

                filled = self.retrying_fill_buf()?;

                if let Some(i) = filled.find_byte(b'\n') {
                    line.extend_from_slice(&filled[..i]);
                    self.consume(i + 1);
                    break line;
                }
            }
        };

        if let Some(b'\r') = line.last() {
            line.pop();
        }
        Ok(line)
    }

    /// Similar to [`BufRead::skip_until`], but continues skipping if the `byte`
    /// does not occur directly after `prefix`.
    ///
    /// It is assumed that `byte` is more rare than `prefix`, and that most
    /// occurrences of `byte` should be preceded by `prefix`.
    ///
    /// For example, this can be used to find `>` at the start of a newline
    /// (after `\n`) for parsing FASTA files.
    #[allow(dead_code)]
    fn skip_until_with_prefix(&mut self, byte: u8, prefix: u8) -> std::io::Result<bool> {
        let mut ends_with_prefix = false;

        loop {
            let filled = self.retrying_fill_buf()?;

            let Some(&last) = filled.last() else {
                // EOF reached without finding pattern
                return Ok(false);
            };

            let used = filled.len();

            if let Some(byte_idx) = filled.find_byte(byte) {
                if let Some(prefix_idx) = byte_idx.checked_sub(1) {
                    if filled[prefix_idx] == prefix {
                        self.consume(byte_idx + 1);
                        return Ok(true);
                    }
                } else if ends_with_prefix {
                    self.consume(1);
                    return Ok(true);
                }

                ends_with_prefix = byte == prefix;
                self.consume(byte_idx + 1);
            } else {
                ends_with_prefix = last == prefix;
                self.consume(used);
            }
        }
    }
}

impl<R: BufRead> BufReadExtension for R {}
