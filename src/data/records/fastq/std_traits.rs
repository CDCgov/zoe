use crate::prelude::*;

impl std::fmt::Display for FastQ {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "@{}\n{}\n+\n{}\n", self.header, self.sequence, self.quality)
    }
}

impl std::fmt::Display for FastQView<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "@{}\n{}\n+\n{}\n", self.header, self.sequence, self.quality)
    }
}

impl std::fmt::Display for FastQViewMut<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "@{}\n{}\n+\n{}\n", self.header, self.sequence, self.quality)
    }
}
