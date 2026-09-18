use std::{fmt::Display, num::ParseIntError, str::FromStr};

/// A wrapper around a `u16` holding the SAM bit flags for alignment records.
///
/// This type by default preserves the unused bits in the `u16`. This impacts
/// equality and hashing. Call [`standardize`] to clear these bits.
///
/// [`standardize`]: Flag::standardize
#[repr(transparent)]
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct Flag(pub u16);

impl Flag {
    /// A [`Flag`] with the unmapped bit set.
    pub(crate) const UNMAPPED: Flag = {
        let mut flag = Flag::new(0);
        flag.set_unmapped();
        flag
    };

    /// Constructs a new [`Flag`] struct by wrapping the given `u16`.
    ///
    /// No standardization of unused bits is performed.
    #[inline]
    #[must_use]
    pub const fn new(flag: u16) -> Self {
        Self(flag)
    }

    /// Clears the unused bits of the `u16` that do not correspond to flags in
    /// the SAM format.
    #[inline]
    pub const fn standardize(&mut self) {
        self.0 &= (1 << 12) - 1;
    }

    /// Returns whether the `0x1` bit of [`Flag`] is set, meaning that the
    /// template has multiple segments in sequencing.
    ///
    /// For example, some programs set this flag when reads are paired.
    #[inline]
    #[must_use]
    pub const fn is_segmented(self) -> bool {
        self.0 & 0x1 > 0
    }

    /// Sets the `0x1` bit of [`Flag`] to indicate that the template has
    /// multiple segments in sequencing.
    ///
    /// For example, some programs set this flag when reads are paired.
    #[inline]
    pub const fn set_segmented(&mut self) {
        self.0 |= 0x1;
    }

    /// Unsets the `0x1` bit of [`Flag`]. When set, this bit indicates that the
    /// template has multiple segments in sequencing.
    #[inline]
    pub const fn unset_segmented(&mut self) {
        self.0 &= !0x1;
    }

    /// Returns whether the `0x2` bit of [`Flag`] is set, meaning that each
    /// segment properly aligned according to the aligner.
    ///
    /// For example, this may be set or unset by a program depending on whether
    /// the orientation and distance between paired reads matches expectations.
    #[inline]
    #[must_use]
    pub const fn is_properly_segmented(self) -> bool {
        self.0 & 0x2 > 0
    }

    /// Sets the `0x2` bit of [`Flag`] to indicate that each segment properly
    /// aligned according to the aligner.
    ///
    /// For example, this may be set or unset by a program depending on whether
    /// the orientation and distance between paired reads matches expectations.
    #[inline]
    pub const fn set_properly_segmented(&mut self) {
        self.0 |= 0x2;
    }

    /// Unsets the `0x2` bit of [`Flag`]. When set, this bit indicates that each
    /// segment properly aligned according to the aligner.
    #[inline]
    pub const fn unset_properly_segmented(&mut self) {
        self.0 &= !0x2;
    }

    /// Returns whether the `0x4` bit of [`Flag`] is set, meaning that the
    /// segment is unmapped.
    ///
    /// For example, this may be set by a program if a read fails to align to a
    /// reference.
    #[inline]
    #[must_use]
    pub const fn is_unmapped(self) -> bool {
        self.0 & 0x4 > 0
    }

    /// Sets the `0x4` bit of [`Flag`] to indicate that the segment is unmapped.
    ///
    /// For example, this may be set by a program if a read fails to align to a
    /// reference.
    #[inline]
    pub const fn set_unmapped(&mut self) {
        self.0 |= 0x4;
    }

    /// Unsets the `0x4` bit of [`Flag`]. When set, this bit indicates that the
    /// segment is unmapped.
    #[inline]
    pub const fn unset_unmapped(&mut self) {
        self.0 &= !0x4;
    }

    /// Returns whether the `0x8` bit of [`Flag`] is set, meaning that the next
    /// segment in the template is unmapped.
    ///
    /// For example, this may be set by a program if the corresponding pair to a
    /// read (in paired end sequencing) is unmapped.
    #[inline]
    #[must_use]
    pub const fn has_unmapped_next_segment(self) -> bool {
        self.0 & 0x8 > 0
    }

    /// Sets the `0x8` bit of [`Flag`] to indicate that the next segment in the
    /// template is unmapped.
    ///
    /// For example, this may be set by a program if the corresponding pair to a
    /// read (in paired end sequencing) is unmapped.
    #[inline]
    pub const fn set_unmapped_next_segment(&mut self) {
        self.0 |= 0x8;
    }

    /// Unsets the `0x8` bit of [`Flag`]. When set, this bit indicates that the
    /// next segment in the template is unmapped.
    #[inline]
    pub const fn unset_unmapped_next_segment(&mut self) {
        self.0 &= !0x8;
    }

    /// Returns whether the `0x10` bit of [`Flag`] is set, meaning that the
    /// `SEQ` field is being reverse complemented.
    ///
    /// For example, for a read aligning to the reverse complement of the
    /// reference, programs typically store the reverse complement of the
    /// sequence so that the alignment corresponds to the forward strand of the
    /// reference. In this case, the program would set this flag.
    #[inline]
    #[must_use]
    pub const fn is_revcomp(self) -> bool {
        self.0 & 0x10 > 0
    }

    /// Sets the `0x10` bit of [`Flag`] to indicate that the `SEQ` field is
    /// being reverse complemented.
    ///
    /// For example, for a read aligning to the reverse complement of the
    /// reference, programs typically store the reverse complement of the
    /// sequence so that the alignment corresponds to the forward strand of the
    /// reference. In this case, the program would set this flag.
    #[inline]
    pub const fn set_revcomp(&mut self) {
        self.0 |= 0x10;
    }

    /// Unsets the `0x10` bit of [`Flag`]. When set, this bit indicates that the
    /// `SEQ` field is being reverse complemented.
    #[inline]
    pub const fn unset_revcomp(&mut self) {
        self.0 &= !0x10;
    }

    /// Returns whether the `0x20` bit of [`Flag`] is set, meaning that the
    /// `SEQ` field of the next segment in the template is being reverse
    /// complemented.
    ///
    /// For example, this may be set by a program if the corresponding pair to a
    /// read (in paired end sequencing) aligned to the reverse complement of the
    /// reference.
    #[inline]
    #[must_use]
    pub const fn has_revcomp_next_segment(self) -> bool {
        self.0 & 0x20 > 0
    }

    /// Sets the `0x20` bit of [`Flag`] to indicate that the `SEQ` field of the
    /// next segment in the template is being reverse complemented.
    ///
    /// For example, this may be set by a program if the corresponding pair to a
    /// read (in paired end sequencing) aligned to the reverse complement of the
    /// reference.
    #[inline]
    pub const fn set_revcomp_next_segment(&mut self) {
        self.0 |= 0x20;
    }

    /// Unsets the `0x20` bit of [`Flag`]. When set, this bit indicates that the
    /// `SEQ` field of the next segment in the template is being reverse
    /// complemented.
    #[inline]
    pub const fn unset_revcomp_next_segment(&mut self) {
        self.0 &= !0x20;
    }

    /// Returns whether the `0x40` bit of [`Flag`] is set, meaning that the
    /// record is the first segment in the template.
    ///
    /// For example, for paired end sequencing, this may be set for the forward
    /// read.
    #[inline]
    #[must_use]
    pub const fn is_first_template(self) -> bool {
        self.0 & 0x40 > 0
    }

    /// Sets the `0x40` bit of [`Flag`] to indicate that the record is the first
    /// segment in the template.
    ///
    /// For example, for paired end sequencing, this may be set for the forward
    /// read.
    #[inline]
    pub const fn set_first_template(&mut self) {
        self.0 |= 0x40;
    }

    /// Unsets the `0x40` bit of [`Flag`]. When set, this bit indicates that the
    /// record is the first segment in the template.
    #[inline]
    pub const fn unset_first_template(&mut self) {
        self.0 &= !0x40;
    }

    /// Returns whether the `0x80` bit of [`Flag`] is set, meaning that the
    /// record is the last segment in the template.
    ///
    /// For example, for paired end sequencing, this may be set for the reverse
    /// read.
    #[inline]
    #[must_use]
    pub const fn is_last_template(self) -> bool {
        self.0 & 0x80 > 0
    }

    /// Sets the `0x80` bit of [`Flag`] to indicate that the record is the last
    /// segment in the template.
    ///
    /// For example, for paired end sequencing, this may be set for the reverse
    /// read.
    #[inline]
    pub const fn set_last_template(&mut self) {
        self.0 |= 0x80;
    }

    /// Unsets the `0x80` bit of [`Flag`]. When set, this bit indicates that the
    /// record is the last segment in the template.
    #[inline]
    pub const fn unset_last_template(&mut self) {
        self.0 &= !0x80;
    }

    /// Returns whether the `0x100` bit of [`Flag`] is set, meaning that the
    /// record is a secondary alignment.
    ///
    /// For example, if a read is chimeric, a program might report both
    /// alignments and mark one as secondary.
    #[inline]
    #[must_use]
    pub const fn is_secondary(self) -> bool {
        self.0 & 0x100 > 0
    }

    /// Sets the `0x100` bit of [`Flag`] to indicate that the record is a
    /// secondary alignment.
    ///
    /// For example, if a read is chimeric, a program might report both
    /// alignments and mark one as secondary.
    #[inline]
    pub const fn set_secondary(&mut self) {
        self.0 |= 0x100;
    }

    /// Unsets the `0x100` bit of [`Flag`]. When set, this bit indicates that
    /// the record is a secondary alignment.
    #[inline]
    pub const fn unset_secondary(&mut self) {
        self.0 &= !0x100;
    }

    /// Returns whether the `0x200` bit of [`Flag`] is set, meaning that the
    /// record does not pass filters, such as platform/vendor quality controls.
    #[inline]
    #[must_use]
    pub const fn fails_qc(self) -> bool {
        self.0 & 0x200 > 0
    }

    /// Sets the `0x200` bit of [`Flag`] to indicate that the record does not
    /// pass filters, such as platform/vendor quality controls.
    pub const fn set_fails_qc(&mut self) {
        self.0 |= 0x200;
    }

    /// Unsets the `0x200` bit of [`Flag`]. When set, this bit indicates that
    /// the record does not pass filters, such as platform/vendor quality
    /// controls.
    pub const fn unset_fails_qc(&mut self) {
        self.0 &= !0x200;
    }

    /// Returns whether the `0x400` bit of [`Flag`] is set, meaning that the
    /// record is a PCR or optical duplicate.
    #[inline]
    #[must_use]
    pub const fn is_duplicate(self) -> bool {
        self.0 & 0x400 > 0
    }

    /// Sets the `0x400` bit of [`Flag`] to indicate that the record is a PCR or
    /// optical duplicate.
    #[inline]
    pub const fn set_duplicate(&mut self) {
        self.0 |= 0x400;
    }

    /// Unsets the `0x400` bit of [`Flag`]. When set, this bit indicates that
    /// the record is a PCR or optical duplicate.
    #[inline]
    pub const fn unset_duplicate(&mut self) {
        self.0 &= !0x400;
    }

    /// Returns whether the `0x800` bit of [`Flag`] is set, meaning that the
    /// record is a supplementary alignment.
    #[inline]
    #[must_use]
    pub const fn is_supplementary(self) -> bool {
        self.0 & 0x800 > 0
    }

    /// Sets the `0x800` bit of [`Flag`] to indicate that the record is a
    /// supplementary alignment.
    #[inline]
    pub const fn set_supplementary(&mut self) {
        self.0 |= 0x800;
    }

    /// Unsets the `0x800` bit of [`Flag`]. When set, this bit indicates that
    /// the record is a supplementary alignment.
    #[inline]
    pub const fn unset_supplementary(&mut self) {
        self.0 &= !0x800;
    }
}

impl Display for Flag {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl FromStr for Flag {
    type Err = ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.parse().map(Flag)
    }
}

#[cfg(test)]
mod test {
    use crate::data::sam::Flag;

    #[test]
    fn segmented() {
        let mut flag = Flag::default();
        assert!(!flag.is_segmented());
        flag.unset_segmented();
        assert_eq!(flag, Flag::default());
        flag.set_segmented();
        assert!(flag.is_segmented());
        assert_eq!(flag.0.count_ones(), 1);
        let mut new_flag = flag;
        new_flag.set_segmented();
        assert_eq!(flag, new_flag);
        flag.unset_segmented();
        assert!(!flag.is_segmented());
        assert_eq!(flag, Flag::default());

        let mut flag = Flag(u16::MAX);
        assert!(flag.is_segmented());
        flag.set_segmented();
        assert_eq!(flag, Flag(u16::MAX));
        flag.unset_segmented();
        assert!(!flag.is_segmented());
        assert_eq!(flag.0.count_ones(), 15);
        let mut new_flag = flag;
        new_flag.unset_segmented();
        assert_eq!(flag, new_flag);
        flag.set_segmented();
        assert!(flag.is_segmented());
        assert_eq!(flag, Flag(u16::MAX));
    }

    #[test]
    fn properly_segmented() {
        let mut flag = Flag::default();
        assert!(!flag.is_properly_segmented());
        flag.unset_properly_segmented();
        assert_eq!(flag, Flag::default());
        flag.set_properly_segmented();
        assert!(flag.is_properly_segmented());
        assert_eq!(flag.0.count_ones(), 1);
        let mut new_flag = flag;
        new_flag.set_properly_segmented();
        assert_eq!(flag, new_flag);
        flag.unset_properly_segmented();
        assert!(!flag.is_properly_segmented());
        assert_eq!(flag, Flag::default());

        let mut flag = Flag(u16::MAX);
        assert!(flag.is_properly_segmented());
        flag.set_properly_segmented();
        assert_eq!(flag, Flag(u16::MAX));
        flag.unset_properly_segmented();
        assert!(!flag.is_properly_segmented());
        assert_eq!(flag.0.count_ones(), 15);
        let mut new_flag = flag;
        new_flag.unset_properly_segmented();
        assert_eq!(flag, new_flag);
        flag.set_properly_segmented();
        assert!(flag.is_properly_segmented());
        assert_eq!(flag, Flag(u16::MAX));
    }

    #[test]
    fn unmapped() {
        let mut flag = Flag::default();
        assert!(!flag.is_unmapped());
        flag.unset_unmapped();
        assert_eq!(flag, Flag::default());
        flag.set_unmapped();
        assert!(flag.is_unmapped());
        assert_eq!(flag.0.count_ones(), 1);
        let mut new_flag = flag;
        new_flag.set_unmapped();
        assert_eq!(flag, new_flag);
        flag.unset_unmapped();
        assert!(!flag.is_unmapped());
        assert_eq!(flag, Flag::default());

        let mut flag = Flag(u16::MAX);
        assert!(flag.is_unmapped());
        flag.set_unmapped();
        assert_eq!(flag, Flag(u16::MAX));
        flag.unset_unmapped();
        assert!(!flag.is_unmapped());
        assert_eq!(flag.0.count_ones(), 15);
        let mut new_flag = flag;
        new_flag.unset_unmapped();
        assert_eq!(flag, new_flag);
        flag.set_unmapped();
        assert!(flag.is_unmapped());
        assert_eq!(flag, Flag(u16::MAX));
    }

    #[test]
    fn unmapped_next_segment() {
        let mut flag = Flag::default();
        assert!(!flag.has_unmapped_next_segment());
        flag.unset_unmapped_next_segment();
        assert_eq!(flag, Flag::default());
        flag.set_unmapped_next_segment();
        assert!(flag.has_unmapped_next_segment());
        assert_eq!(flag.0.count_ones(), 1);
        let mut new_flag = flag;
        new_flag.set_unmapped_next_segment();
        assert_eq!(flag, new_flag);
        flag.unset_unmapped_next_segment();
        assert!(!flag.has_unmapped_next_segment());
        assert_eq!(flag, Flag::default());

        let mut flag = Flag(u16::MAX);
        assert!(flag.has_unmapped_next_segment());
        flag.set_unmapped_next_segment();
        assert_eq!(flag, Flag(u16::MAX));
        flag.unset_unmapped_next_segment();
        assert!(!flag.has_unmapped_next_segment());
        assert_eq!(flag.0.count_ones(), 15);
        let mut new_flag = flag;
        new_flag.unset_unmapped_next_segment();
        assert_eq!(flag, new_flag);
        flag.set_unmapped_next_segment();
        assert!(flag.has_unmapped_next_segment());
        assert_eq!(flag, Flag(u16::MAX));
    }

    #[test]
    fn revcomp() {
        let mut flag = Flag::default();
        assert!(!flag.is_revcomp());
        flag.unset_revcomp();
        assert_eq!(flag, Flag::default());
        flag.set_revcomp();
        assert!(flag.is_revcomp());
        assert_eq!(flag.0.count_ones(), 1);
        let mut new_flag = flag;
        new_flag.set_revcomp();
        assert_eq!(flag, new_flag);
        flag.unset_revcomp();
        assert!(!flag.is_revcomp());
        assert_eq!(flag, Flag::default());

        let mut flag = Flag(u16::MAX);
        assert!(flag.is_revcomp());
        flag.set_revcomp();
        assert_eq!(flag, Flag(u16::MAX));
        flag.unset_revcomp();
        assert!(!flag.is_revcomp());
        assert_eq!(flag.0.count_ones(), 15);
        let mut new_flag = flag;
        new_flag.unset_revcomp();
        assert_eq!(flag, new_flag);
        flag.set_revcomp();
        assert!(flag.is_revcomp());
        assert_eq!(flag, Flag(u16::MAX));
    }

    #[test]
    fn revcomp_next_segment() {
        let mut flag = Flag::default();
        assert!(!flag.has_revcomp_next_segment());
        flag.unset_revcomp_next_segment();
        assert_eq!(flag, Flag::default());
        flag.set_revcomp_next_segment();
        assert!(flag.has_revcomp_next_segment());
        assert_eq!(flag.0.count_ones(), 1);
        let mut new_flag = flag;
        new_flag.set_revcomp_next_segment();
        assert_eq!(flag, new_flag);
        flag.unset_revcomp_next_segment();
        assert!(!flag.has_revcomp_next_segment());
        assert_eq!(flag, Flag::default());

        let mut flag = Flag(u16::MAX);
        assert!(flag.has_revcomp_next_segment());
        flag.set_revcomp_next_segment();
        assert_eq!(flag, Flag(u16::MAX));
        flag.unset_revcomp_next_segment();
        assert!(!flag.has_revcomp_next_segment());
        assert_eq!(flag.0.count_ones(), 15);
        let mut new_flag = flag;
        new_flag.unset_revcomp_next_segment();
        assert_eq!(flag, new_flag);
        flag.set_revcomp_next_segment();
        assert!(flag.has_revcomp_next_segment());
        assert_eq!(flag, Flag(u16::MAX));
    }

    #[test]
    fn first_template() {
        let mut flag = Flag::default();
        assert!(!flag.is_first_template());
        flag.unset_first_template();
        assert_eq!(flag, Flag::default());
        flag.set_first_template();
        assert!(flag.is_first_template());
        assert_eq!(flag.0.count_ones(), 1);
        let mut new_flag = flag;
        new_flag.set_first_template();
        assert_eq!(flag, new_flag);
        flag.unset_first_template();
        assert!(!flag.is_first_template());
        assert_eq!(flag, Flag::default());

        let mut flag = Flag(u16::MAX);
        assert!(flag.is_first_template());
        flag.set_first_template();
        assert_eq!(flag, Flag(u16::MAX));
        flag.unset_first_template();
        assert!(!flag.is_first_template());
        assert_eq!(flag.0.count_ones(), 15);
        let mut new_flag = flag;
        new_flag.unset_first_template();
        assert_eq!(flag, new_flag);
        flag.set_first_template();
        assert!(flag.is_first_template());
        assert_eq!(flag, Flag(u16::MAX));
    }

    #[test]
    fn last_template() {
        let mut flag = Flag::default();
        assert!(!flag.is_last_template());
        flag.unset_last_template();
        assert_eq!(flag, Flag::default());
        flag.set_last_template();
        assert!(flag.is_last_template());
        assert_eq!(flag.0.count_ones(), 1);
        let mut new_flag = flag;
        new_flag.set_last_template();
        assert_eq!(flag, new_flag);
        flag.unset_last_template();
        assert!(!flag.is_last_template());
        assert_eq!(flag, Flag::default());

        let mut flag = Flag(u16::MAX);
        assert!(flag.is_last_template());
        flag.set_last_template();
        assert_eq!(flag, Flag(u16::MAX));
        flag.unset_last_template();
        assert!(!flag.is_last_template());
        assert_eq!(flag.0.count_ones(), 15);
        let mut new_flag = flag;
        new_flag.unset_last_template();
        assert_eq!(flag, new_flag);
        flag.set_last_template();
        assert!(flag.is_last_template());
        assert_eq!(flag, Flag(u16::MAX));
    }

    #[test]
    fn secondary() {
        let mut flag = Flag::default();
        assert!(!flag.is_secondary());
        flag.unset_secondary();
        assert_eq!(flag, Flag::default());
        flag.set_secondary();
        assert!(flag.is_secondary());
        assert_eq!(flag.0.count_ones(), 1);
        let mut new_flag = flag;
        new_flag.set_secondary();
        assert_eq!(flag, new_flag);
        flag.unset_secondary();
        assert!(!flag.is_secondary());
        assert_eq!(flag, Flag::default());

        let mut flag = Flag(u16::MAX);
        assert!(flag.is_secondary());
        flag.set_secondary();
        assert_eq!(flag, Flag(u16::MAX));
        flag.unset_secondary();
        assert!(!flag.is_secondary());
        assert_eq!(flag.0.count_ones(), 15);
        let mut new_flag = flag;
        new_flag.unset_secondary();
        assert_eq!(flag, new_flag);
        flag.set_secondary();
        assert!(flag.is_secondary());
        assert_eq!(flag, Flag(u16::MAX));
    }

    #[test]
    fn fails_qc() {
        let mut flag = Flag::default();
        assert!(!flag.fails_qc());
        flag.unset_fails_qc();
        assert_eq!(flag, Flag::default());
        flag.set_fails_qc();
        assert!(flag.fails_qc());
        assert_eq!(flag.0.count_ones(), 1);
        let mut new_flag = flag;
        new_flag.set_fails_qc();
        assert_eq!(flag, new_flag);
        flag.unset_fails_qc();
        assert!(!flag.fails_qc());
        assert_eq!(flag, Flag::default());

        let mut flag = Flag(u16::MAX);
        assert!(flag.fails_qc());
        flag.set_fails_qc();
        assert_eq!(flag, Flag(u16::MAX));
        flag.unset_fails_qc();
        assert!(!flag.fails_qc());
        assert_eq!(flag.0.count_ones(), 15);
        let mut new_flag = flag;
        new_flag.unset_fails_qc();
        assert_eq!(flag, new_flag);
        flag.set_fails_qc();
        assert!(flag.fails_qc());
        assert_eq!(flag, Flag(u16::MAX));
    }

    #[test]
    fn duplicate() {
        let mut flag = Flag::default();
        assert!(!flag.is_duplicate());
        flag.unset_duplicate();
        assert_eq!(flag, Flag::default());
        flag.set_duplicate();
        assert!(flag.is_duplicate());
        assert_eq!(flag.0.count_ones(), 1);
        let mut new_flag = flag;
        new_flag.set_duplicate();
        assert_eq!(flag, new_flag);
        flag.unset_duplicate();
        assert!(!flag.is_duplicate());
        assert_eq!(flag, Flag::default());

        let mut flag = Flag(u16::MAX);
        assert!(flag.is_duplicate());
        flag.set_duplicate();
        assert_eq!(flag, Flag(u16::MAX));
        flag.unset_duplicate();
        assert!(!flag.is_duplicate());
        assert_eq!(flag.0.count_ones(), 15);
        let mut new_flag = flag;
        new_flag.unset_duplicate();
        assert_eq!(flag, new_flag);
        flag.set_duplicate();
        assert!(flag.is_duplicate());
        assert_eq!(flag, Flag(u16::MAX));
    }

    #[test]
    fn supplementary() {
        let mut flag = Flag::default();
        assert!(!flag.is_supplementary());
        flag.unset_supplementary();
        assert_eq!(flag, Flag::default());
        flag.set_supplementary();
        assert!(flag.is_supplementary());
        assert_eq!(flag.0.count_ones(), 1);
        let mut new_flag = flag;
        new_flag.set_supplementary();
        assert_eq!(flag, new_flag);
        flag.unset_supplementary();
        assert!(!flag.is_supplementary());
        assert_eq!(flag, Flag::default());

        let mut flag = Flag(u16::MAX);
        assert!(flag.is_supplementary());
        flag.set_supplementary();
        assert_eq!(flag, Flag(u16::MAX));
        flag.unset_supplementary();
        assert!(!flag.is_supplementary());
        assert_eq!(flag.0.count_ones(), 15);
        let mut new_flag = flag;
        new_flag.unset_supplementary();
        assert_eq!(flag, new_flag);
        flag.set_supplementary();
        assert!(flag.is_supplementary());
        assert_eq!(flag, Flag(u16::MAX));
    }
}
