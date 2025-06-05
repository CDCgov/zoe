pub(crate) trait ChopLineBreak {
    /// Removes a trailing `\n` or `\r\n`
    fn chop_line_break(&mut self);
}

impl ChopLineBreak for Vec<u8> {
    #[inline]
    fn chop_line_break(&mut self) {
        if self.ends_with(b"\n") {
            self.pop();

            if self.ends_with(b"\r") {
                self.pop();
            }
        }
    }
}

pub(crate) trait StripLineBreak {
    /// Returns a new slice without a trailing `\n` or `\r\n`
    fn strip_line_break(&self) -> Self;
}

impl<'a> StripLineBreak for &'a [u8] {
    #[inline]
    fn strip_line_break(&self) -> &'a [u8] {
        let mut out = *self;
        if out.ends_with(b"\n") {
            out = &out[..out.len() - 1];

            if out.ends_with(b"\r") {
                out = &out[..out.len() - 1];
            }
        }
        out
    }
}
