use crate::search::ByteSubstring;
use std::io::{BufRead, ErrorKind};

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
}

impl<R: BufRead> BufReadExtension for R {}
