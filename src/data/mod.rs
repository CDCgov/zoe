//! ## Data import, export, and manipulation functions.
//!
//! ## IUPAC Standards
//!
//! For nucleotides and proteins, *Zoe* implements [IUPAC
//! definitions](https://www.bioinformatics.org/sms/iupac.html) that include
//! ambiguous base calls.
//!
//! For amino acid scoring matrices, the following ambiguous amino acids are
//! also included:
//! - `B`, which can be either `D` or `N` (Aspartic acid or Asparagine
//!   respectively)
//! - `Z`, which can be either `E` or `Q` (Glutamic acid or Glutamine
//!   respectively)
//! - `J`, which can be either `I` or `L` (Isoleucine or Leucine respectively)
//! - `X`, which represents an unknown position
//!
//! `B` and `Z` are common due to the pairs sharing the same carbon backbone,
//! and since a deamidation reaction may occur during sequencing which cause
//! ambiguity (see [here](https://biology.stackexchange.com/a/41350)). `J` is
//! common because the pair are constitutional isomers and have similar mass
//! spectrums, which can make them difficult to distinguish in some cases (e.g.,
//! [this paper](https://pubs.acs.org/doi/10.1021/ac020422m#)).
//!
//! ## Views
//!
//! Many of the data type provided by *Zoe* have versions holding owned data as
//! well as versions holding references. The latter are called *views*, and are
//! useful when performing operations on a subsequence of the original data. For
//! example, [`Nucleotides`] (which is a wrapper around [`Vec<u8>`]) has the
//! corresponding types [`NucleotidesView`] and [`NucleotidesViewMut`], which
//! are wrappers around immutable and mutable byte slices respectively.
//!
//! A view can be constructed directly from an existing slice, such as
//! [`NucleotidesView::from_bytes_unchecked`], but more often a view is created
//! from an owned instance of the data by calling [`as_view`] or
//! [`as_view_mut`]. If only a range of the data is desired, then [`slice`] and
//! [`slice_mut`] are used.
//!
//! A view can be directly displayed/printed, or it can be copied into owned
//! data again by using [`to_owned_data`].
//!
//! For example:
//! ```
//! # use zoe::prelude::*;
//! let owned_sequence: Nucleotides = b"GGCCACCAAGGCCA".into();
//!
//! let view_of_sequence = owned_sequence.as_view();
//! assert_eq!(view_of_sequence.as_bytes(), b"GGCCACCAAGGCCA");
//! assert_eq!(view_of_sequence.gc_content(), 10);
//!
//! let slice_of_middle = owned_sequence.slice(4..9);
//! assert_eq!(slice_of_middle.as_bytes(), b"ACCAA");
//! assert_eq!(slice_of_middle.gc_content(), 2);
//!
//! let slice_of_middle_again = slice_of_middle.slice(2..);
//! assert_eq!(slice_of_middle_again.as_bytes(), b"CAA");
//!
//! let view_to_owned = slice_of_middle.to_owned();
//! assert_eq!(view_to_owned, b"ACCAA".into());
//!
//! let mut owned_sequence: Nucleotides = b"GGCCACCAAGGCCA".into();
//! let mut mutable_slice = owned_sequence.slice_mut(1..4);
//! mutable_slice.as_mut_bytes().fill(b'T');
//! assert_eq!(owned_sequence.as_bytes(), b"GTTTACCAAGGCCA");
//! ```
//!
//! Views can also be re-sliced in-place using the [`restrict`] method. This is
//! a useful way to avoid needing to have extra let bindings. For example:
//! ```
//! # use zoe::prelude::*;
//! let mut owned_sequence: Nucleotides = b"GGCCACCAAGGCCA".into();
//! let mut mutable_slice = owned_sequence.as_view_mut();
//! mutable_slice.restrict(5..);
//! mutable_slice.restrict(..2);
//! assert_eq!(mutable_slice.as_bytes(), b"CC");
//! ```
//!
//! One exception to the above is [`Cigar`] strings. [`CigarView`] and
//! [`CigarViewMut`] are available and can be created with [`as_view`] and
//! [`as_view_mut`]. However, slicing and restricting are not possible for these
//! types, since they are not sequence data.
//!
//! ## IO Errors in *Zoe*
//!
//! As a library, *Zoe* aims to avoid making assumptions on the style of error
//! handling chosen by users, in particular by not adopting any error handling
//! crate as a dependency.
//!
//! For specific applications, *Zoe* has enum-style error types such as
//! [`QueryProfileError`] or [`KmerError`], which the user can match on or
//! display. For working with files and record types, however, *Zoe* elects to
//! use [`std::io::Error`], allowing for system IO errors to be propagated and
//! function-specific error messages to be represented with
//! [`ErrorKind::InvalidData`] or [`ErrorKind::Other`].
//!
//! IO failures are assumed to be rare by *Zoe*, and hence the crate will
//! automatically add the file path to the error messages (such as
//! [`FastQReader::from_filename`]). If in doubt, check the `Errors` section of
//! a function's documentation to determine what information is automatically
//! added.
//!
//! When *Zoe* adds context to an error message, it will store the original
//! error message so that it is accessible using [`Error::source`]. Many error
//! handling crates such as `anyhow` automatically will display this information
//! in the backtrace. Or, *Zoe* offers [`unwrap_or_fail`] and [`unwrap_or_die`]
//! which will include this information when unwrapping an error.
//!
//! [`Nucleotides`]: types::nucleotides::Nucleotides
//! [`NucleotidesView`]: types::nucleotides::NucleotidesView
//! [`NucleotidesViewMut`]: types::nucleotides::NucleotidesViewMut
//! [`NucleotidesView::from_bytes_unchecked`]:
//!     types::nucleotides::NucleotidesView::from_bytes_unchecked
//! [`Cigar`]: types::cigar::Cigar
//! [`CigarView`]: types::cigar::CigarView
//! [`CigarViewMut`]: types::cigar::CigarViewMut
//! [`as_view`]: crate::prelude::DataOwned::as_view
//! [`as_view_mut`]: crate::prelude::DataOwned::as_view_mut
//! [`slice`]: crate::prelude::Slice::slice
//! [`slice_mut`]: crate::prelude::SliceMut::slice_mut
//! [`to_owned_data`]: crate::prelude::DataView::to_owned_data
//! [`restrict`]: crate::prelude::Restrict::restrict
//! [`QueryProfileError`]: crate::data::err::QueryProfileError
//! [`KmerError`]: crate::kmer::KmerError
//! [`ErrorKind::InvalidData`]: std::io::ErrorKind::InvalidData
//! [`ErrorKind::Other`]: std::io::ErrorKind::Other
//! [`FastQReader::from_filename`]: crate::data::records::fastq::FastQReader
//! [`Error::source`]: std::error::Error::source
//! [`unwrap_or_fail`]: crate::data::err::OrFail::unwrap_or_fail
//! [`unwrap_or_die`]: crate::data::err::OrFail::unwrap_or_die

#[cfg(feature = "fuzzing")]
pub mod arbitrary;
/// A module with error types and convenience traits for handling [`Result`].
pub mod err;
pub mod matrices;
/// A module for records types--usually for I/O--that are structures of other
/// more primitive types.
pub mod records;
/// A module for storing more fundamental types, like
/// [`Nucleotides`](self::types::nucleotides::Nucleotides) and
/// [`AminoAcids`][self::types::amino_acids::AminoAcids].
pub mod types;
/// A module containing functionality needed for working with views.
pub mod views;

/// A private module for helper alphabets, maps, and matrices that can be used
/// within public methods.
pub(crate) mod constants;
/// A private module for helper extension traits.
pub(crate) mod extension;

/// Used for type validation
mod validation;
/// Contains the macros [`define_whichever`] and [`impl_traits`]
///
/// [`define_whichever`]: crate::define_whichever
/// [`impl_traits`]: crate::impl_traits
mod whichever;

pub use constants::mappings::{
    AA_ALL_AMBIG_PROFILE_MAP_WITH_STOP, AA_PROFILE_MAP, AA_UNAMBIG_PROFILE_MAP, ByteIndexMap, DNA_PROFILE_MAP,
    DNA_UNAMBIG_PROFILE_MAP, StdGeneticCode,
};
pub use matrices::WeightMatrix;
pub use records::{fasta, fastq, sam};
pub use types::{amino_acids, cigar, nucleotides, phred};
pub use validation::{CheckSequence, Recode, RetainSequence, StdForSequences};

pub(crate) use constants::{alphas, mappings};
pub(crate) use extension::{array_types, byte_types, id_types, vec_types};

#[cfg(test)]
mod tests {
    use crate::{
        data::cigar::Cigar,
        prelude::{
            AminoAcids, AminoAcidsView, AminoAcidsViewMut, Nucleotides, NucleotidesView, NucleotidesViewMut, QualityScores,
            QualityScoresView, QualityScoresViewMut,
        },
    };

    #[test]
    fn test_from() {
        const SEQ: &[u8; 14] = b"10M3X5D10M3X5D";
        let mut mut_arr1 = *SEQ;
        let mut mut_arr2 = *SEQ;

        macro_rules! test_from_impls_view_mut {
            ($t:ty) => {
                let mut_byte_slice = mut_arr1.as_mut_slice();
                let arr_mut_ref = &mut mut_arr2;
                let mut owned = *SEQ;
                let expected = <$t>::from(&mut owned);
                assert_eq!(<$t>::from(mut_byte_slice), expected);
                assert_eq!(<$t>::from(arr_mut_ref), expected);
            };
        }

        macro_rules! test_from_impls_view {
            ($t:ty) => {
                let byte_slice = SEQ.as_slice();
                let arr_ref = SEQ;
                let mut owned = *SEQ;
                let expected = <$t>::from(&mut owned);
                assert_eq!(<$t>::from(byte_slice), expected);
                assert_eq!(<$t>::from(arr_ref), expected);
                test_from_impls_view_mut!($t);
            };
        }

        macro_rules! test_from_impls_owned {
            ($t:ty) => {
                let string = String::from_utf8_lossy(SEQ.as_slice()).to_string();
                let vec = SEQ.to_vec();
                let arr = *SEQ;
                let mut owned = *SEQ;
                let expected = <$t>::from(&mut owned);
                assert_eq!(<$t>::from(string), expected);
                assert_eq!(<$t>::from(vec), expected);
                assert_eq!(<$t>::from(arr), expected);
                test_from_impls_view!($t);
            };
        }

        macro_rules! test_try_from_impls_view_mut {
            ($t:ty) => {
                let mut_byte_slice = mut_arr1.as_mut_slice();
                let arr_mut_ref = &mut mut_arr2;
                let mut owned = *SEQ;
                let expected = <$t>::try_from(&mut owned).unwrap();
                assert_eq!(<$t>::try_from(mut_byte_slice).unwrap(), expected);
                assert_eq!(<$t>::try_from(arr_mut_ref).unwrap(), expected);
            };
        }

        macro_rules! test_try_from_impls_view {
            ($t:ty) => {
                let byte_slice = SEQ.as_slice();
                let arr_ref = SEQ;
                let mut owned = *SEQ;
                let expected = <$t>::try_from(&mut owned).unwrap();
                assert_eq!(<$t>::try_from(byte_slice).unwrap(), expected);
                assert_eq!(<$t>::try_from(arr_ref).unwrap(), expected);
                test_try_from_impls_view_mut!($t);
            };
        }

        macro_rules! test_try_from_impls_owned {
            ($t:ty) => {
                let string = String::from_utf8_lossy(SEQ.as_slice()).to_string();
                let vec = SEQ.to_vec();
                let arr = *SEQ;
                let mut owned = *SEQ;
                let expected = <$t>::try_from(&mut owned).unwrap();
                assert_eq!(<$t>::try_from(string).unwrap(), expected);
                assert_eq!(<$t>::try_from(vec).unwrap(), expected);
                assert_eq!(<$t>::try_from(arr).unwrap(), expected);
                test_try_from_impls_view!($t);
            };
        }

        test_from_impls_owned!(Nucleotides);
        test_from_impls_view!(NucleotidesView);
        test_from_impls_view_mut!(NucleotidesViewMut);
        test_from_impls_owned!(AminoAcids);
        test_from_impls_view!(AminoAcidsView);
        test_from_impls_view_mut!(AminoAcidsViewMut);
        test_try_from_impls_owned!(QualityScores);
        test_try_from_impls_view!(QualityScoresView);
        test_try_from_impls_view_mut!(QualityScoresViewMut);
        test_try_from_impls_owned!(Cigar);
    }
}
