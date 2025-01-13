pub(crate) trait ChopLineBreak {
    /// Removes trailing `\n` or `\r\n`
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
