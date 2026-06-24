use crate::data::{
    ByteValidator,
    arbitrary::{ArbitrarySpecs, ByteSet, ByteSpecs, Case, StringSpecs},
    fasta::generic::{Fasta, FastaAnnot},
};
use arbitrary::{Arbitrary, Result, Unstructured};

impl<'a, S> Arbitrary<'a> for Fasta<S>
where
    S: Arbitrary<'a>,
{
    /// Generates an arbitrary [`Fasta`] record from the given unstructured
    /// data.
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let header = String::arbitrary(u)?.replace('>', " ");
        let sequence = S::arbitrary(u)?;

        Ok(Fasta { header, sequence })
    }
}

impl<'a, M, S> Arbitrary<'a> for FastaAnnot<M, S>
where
    M: Arbitrary<'a>,
    S: Arbitrary<'a>,
{
    /// Generates an arbitrary [`FastaAnnot`] record from the given unstructured
    /// data.
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let header = String::arbitrary(u)?.replace('>', " ");
        let sequence = S::arbitrary(u)?;
        let annot = M::arbitrary(u)?;

        Ok(FastaAnnot { header, sequence, annot })
    }
}

/// Specifications for generating an arbitrary [`Fasta`] record.
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct FastaSpecs<S> {
    /// The specifications for generating the header.
    pub header_specs: StringSpecs,

    /// The specifications for generating the sequence.
    pub sequence_specs: S,
}

impl<S> FastaSpecs<S>
where
    S: From<ByteSpecs>,
{
    /// A constructor for [`FastaSpecs`] that ensures the generated [`Fasta`]
    /// can be parsed by [`FastaReader`].
    ///
    /// This ensures that the header and sequence do not contain `\n` or `>`.
    ///
    /// [`FastaReader`]: crate::prelude::FastaReader
    #[inline]
    #[must_use]
    pub fn parsable_fasta() -> Self {
        const BYTE_SET: ByteSet = ByteSet::Custom(
            &ByteValidator::none()
                .add_range(0..=127)
                .remove(b"\n>")
                .generate_alphabet::<126>(),
        );

        Self {
            header_specs:   StringSpecs {
                set:  BYTE_SET,
                case: Case::Any,
            },
            sequence_specs: ByteSpecs {
                set:  BYTE_SET,
                case: Case::Any,
            }
            .into(),
        }
    }
}

impl<'a, S> ArbitrarySpecs<'a> for FastaSpecs<S>
where
    S: ArbitrarySpecs<'a>,
{
    type Output = Fasta<S::Output>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let header = self.header_specs.make_arbitrary(u)?;
        let sequence = self.sequence_specs.make_arbitrary(u)?;

        Ok(Fasta { header, sequence })
    }
}

/// Specifications for generating an arbitrary [`FastaAnnot`] record.
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct FastaAnnotSpecs<M, S> {
    /// The specifications for generating the header.
    pub header_specs: StringSpecs,

    /// The specifications for generating the sequence.
    pub sequence_specs: S,

    /// The specifications for generating the metadata.
    pub metadata_specs: M,
}

impl<'a, M, S> ArbitrarySpecs<'a> for FastaAnnotSpecs<M, S>
where
    S: ArbitrarySpecs<'a>,
    M: ArbitrarySpecs<'a>,
{
    type Output = FastaAnnot<M::Output, S::Output>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let header = self.header_specs.make_arbitrary(u)?;
        let sequence = self.sequence_specs.make_arbitrary(u)?;
        let annot = self.metadata_specs.make_arbitrary(u)?;

        Ok(FastaAnnot { header, sequence, annot })
    }
}
