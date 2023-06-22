pub(crate) trait ByteAvg {
    fn avg(self, other: u8) -> u8;
}

#[allow(clippy::cast_possible_truncation)]
impl ByteAvg for u8 {
    fn avg(self, other: u8) -> u8 {
        ((u16::from(self) + u16::from(other)) / 2u16) as u8
    }
}

#[allow(clippy::wrong_self_convention)]
pub(crate) trait IsBase {
    fn is_base(self) -> bool;
    fn is_not_base(self) -> bool;
}

impl IsBase for u8 {
    #[inline]
    fn is_base(self) -> bool {
        self == b'A' || self == b'G' || self == b'T' || self == b'C'
    }

    #[inline]
    fn is_not_base(self) -> bool {
        !self.is_base()
    }
}
