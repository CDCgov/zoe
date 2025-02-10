//! ## Data import, export, and manipulation functions.
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
//! [`Nucleotides`]: types::nucleotides::Nucleotides
//! [`NucleotidesView`]: types::nucleotides::NucleotidesView
//! [`NucleotidesViewMut`]: types::nucleotides::NucleotidesViewMut
//! [`NucleotidesView::from_bytes_unchecked`]:
//!     types::nucleotides::NucleotidesView::from_bytes_unchecked
//! [`as_view`]: crate::prelude::DataOwned::as_view
//! [`as_view_mut`]: crate::prelude::DataOwned::as_view_mut
//! [`slice`]: crate::prelude::Slice::slice
//! [`slice_mut`]: crate::prelude::SliceMut::slice_mut
//! [`to_owned_data`]: crate::prelude::DataView::to_owned_data
//! [`restrict`]: crate::prelude::DataView::restrict

/// A module with error types and convenience traits for handling [`Result`].
pub mod err;
/// A module for records types--usually for I/O--that are structures of other
/// more primitive types.
pub mod records;
/// A module for storing more fundamental types, like
/// [`Nucleotides`](self::types::nucleotides::Nucleotides) and
/// [`AminoAcids`][self::types::amino_acids::AminoAcids].
pub mod types;
/// A module containing the traits needed for working with views.
pub mod view_traits;

/// A private module for helper alphabets, maps, and matrices that can be used
/// within public methods.
pub(crate) mod constants;
/// A private module for helper extension traits.
pub(crate) mod extension;

/// Used for type validation
mod validation;

pub use constants::{
    mappings::{ByteIndexMap, DNA_PROFILE_MAP, StdGeneticCode},
    matrices::WeightMatrix,
};
pub use records::{fasta, fastq, sam};
pub use types::{amino_acids, cigar, nucleotides, phred};
pub use validation::{CheckSequence, Recode, RetainSequence, StdForSequences};

pub(crate) use constants::{alphas, mappings, matrices};
pub(crate) use extension::{array_types, byte_types, id_types, vec_types};
