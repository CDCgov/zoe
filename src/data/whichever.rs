/// A macro to define enums whose variants all implement particular traits.
///
/// Specifically, this macro is able to define enums whose variants each contain
/// a single type, where the types all implement [`Read`], [`Write`], or
/// [`Iterator`]. This macro is then able to automatically implement that trait
/// for the enum as a whole.
///
/// This is similar to `Either`, but with any number of variants and where the
/// types are known (not generic).
///
/// An invocation will contain:
/// 1. The enum definition, with any number of outer attributes, an optional
///    visibility specifier, and the variants
/// 2. The traits to implement, using `impl Trait {}`. Currently, we support
///    [`Read`], [`Write`], and [`Iterator`].
///
/// ## Examples
///
/// For implementing [`Read`]:
/// ```
/// use std::io::Read;
/// use zoe::define_whichever;
///
/// define_whichever!{
///     #[doc = "An enum representing the allowable input types"]
///     pub enum AnyInput {
///         File(std::fs::File),
///         Stdin(std::io::Stdin),
///         Pipe(std::io::PipeReader),
///     }
///
///     impl Read for AnyInput {}
/// }
/// ```
///
/// For implementing [`Iterator`]:
/// ```
/// #![feature(try_trait_v2)]
/// use zoe::define_whichever;
///
/// define_whichever!{
///     #[doc = "An enum representing iterators yielding the same value several times"]
///     pub enum AnyIterator<T> {
///         Once(std::iter::Once<T>),
///         RepeatN(std::iter::RepeatN<T>),
///         Repeat(std::iter::Repeat<T>),
///     }
///
///     impl<T: Clone> Iterator for AnyIterator<T> {
///         type Item = T;
///     }
/// }
/// ```
///
/// ## Notes
///
/// When implementing [`Iterator`], the generics `F`, `B`, `R`, `P`, `S`, and
/// `K` are reserved. Also, the `try_trait_v2` feature must be enabled.
///
/// ## Acknowledgements
///
/// The syntax for matching generics was developed by
/// [pin-project-lite](https://github.com/taiki-e/pin-project-lite/tree/main).
///
/// [`Read`]: std::io::Read
/// [`Write`]: std::io::Write
#[macro_export]
macro_rules! define_whichever {
    (
        $(#[$meta:meta])*
        $vis:vis enum $struct_name:ident $(<
            $( $lifetime:lifetime $(: $lifetime_bound:lifetime)? ),* $(,)?
            $( $generics:ident
                $(: $generics_bound:path)?
                $(: ?$generics_unsized_bound:path)?
                $(: $generics_lifetime_bound:lifetime)?
            ),*
        >)? {
            $(
                $(#[$variant_meta:meta])*
                $variant:ident($ty:ty)
            ),+
            $(,)?
        }

        $($impl_stuff:tt)*
    ) => {
        macro_rules! match_macro {
            ($value:expr, $pattern:pat => $result:expr) => {
                match $value {
                    $(
                        $struct_name::$variant($pattern) => $result,
                    )+
                }
            };
        }

        $(#[$meta])*
        $vis enum $struct_name $(<
            $( $lifetime $(: $lifetime_bound)? ),*
            $( $generics
                $(: $generics_bound)?
                $(: ?$generics_unsized_bound)?
                $(: $generics_lifetime_bound)?
            ),*
        >)? {
            $(
                $(#[$variant_meta])*
                $variant($ty),
            )+
        }

        $crate::impl_traits!{@match_stmt $($impl_stuff)*}

    };
}

/// A macro to aid in implementing traits for wrapper types. This macro is also
/// used internally by [`define_whichever`].
///
/// This macro currently supports implementing [`Read`], [`Write`], and
/// [`Iterator`].
///
/// ## Examples
///
/// For implementing [`Read`]:
/// ```
/// use std::io::Read;
/// use zoe::impl_traits;
///
/// pub struct MyReader(std::io::BufReader<std::fs::File>);
///
/// impl_traits!{
///     impl Read for MyReader {}
/// }
/// ```
///
/// For implementing [`Iterator`]:
/// ```
/// #![feature(try_trait_v2)]
/// use zoe::impl_traits;
///
/// pub struct MyIterator<T, M>(std::iter::Map<std::iter::RepeatN<T>, M>);
///
/// impl_traits!{
///     impl<T: Clone, Q, M: FnMut(T) -> Q> Iterator for MyIterator<T, M> {
///         type Item = Q;
///     }
/// }
/// ```
///
/// ## Notes
///
/// When implementing [`Iterator`], the generics `F`, `B`, `R`, `P`, `S`, and
/// `K` are reserved. Also, the `try_trait_v2` feature must be enabled.
///
/// ## Acknowledgements
///
/// The syntax for matching generics was developed by
/// [pin-project-lite](https://github.com/taiki-e/pin-project-lite/tree/main).
///
/// [`Read`]: std::io::Read
/// [`Write`]: std::io::Write
#[macro_export]
macro_rules! impl_traits {
    (impl $($tt:tt)*) => {
        $crate::impl_traits!{@wrapper_type impl $($tt)*}
    };

    (@$dispatch:tt
        $(
            impl $(<
                $( $lifetime:lifetime $(: $lifetime_bound:lifetime)? ),* $(,)?
                $( $generics:ident
                    $(: $generics_bound:path)?
                    $(: ?$generics_unsized_bound:path)?
                    $(: $generics_lifetime_bound:lifetime)?
                ),*
            >)? $trait:ident for $self_ty:ty
            $(where
                $( $where_clause_ty:ty
                    $(: $where_clause_bound:path)?
                    $(: ?$where_clause_unsized_bound:path)?
                    $(: $where_clause_lifetime_bound:lifetime)?
                ),* $(,)?
            )?
            {
                $($tt:tt)*
            }
        )*
    ) => {
        $(
            impl $(<
                $( $lifetime $(: $lifetime_bound)? ),*
                $( $generics
                    $(: $generics_bound)?
                    $(: ?$generics_unsized_bound)?
                    $(: $generics_lifetime_bound)?
                ),*
            >)? $trait for $self_ty
            $(where
                $( $where_clause_ty
                    $(: $where_clause_bound)?
                    $(: ?$where_clause_unsized_bound)?
                    $(: $where_clause_lifetime_bound)?
                ),* $(,)?
            )?
            {
                $($tt)*
                $crate::impl_traits!{@$dispatch $trait}
            }
        )*
    };

    (@$dispatch:tt Read) => {
        #[inline]
        fn read(&mut self, buf: &mut [u8]) -> ::std::io::Result<usize> {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.read(buf))
        }

        #[inline]
        fn read_vectored(&mut self, bufs: &mut [::std::io::IoSliceMut<'_>]) -> ::std::io::Result<usize> {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.read_vectored(bufs))
        }

        #[inline]
        fn read_to_end(&mut self, buf: &mut ::std::vec::Vec<u8>) -> ::std::io::Result<usize> {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.read_to_end(buf))
        }

        #[inline]
        fn read_to_string(&mut self, buf: &mut ::std::string::String) -> ::std::io::Result<usize> {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.read_to_string(buf))
        }

        #[inline]
        fn read_exact(&mut self, buf: &mut [u8]) -> ::std::io::Result<()> {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.read_exact(buf))
        }
    };

    (@$dispatch:tt Write) => {
        #[inline]
        fn write(&mut self, buf: &[u8]) -> ::std::io::Result<usize> {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.write(buf))
        }

        #[inline]
        fn flush(&mut self) -> ::std::io::Result<()> {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.flush())
        }

        #[inline]
        fn write_vectored(&mut self, bufs: &[::std::io::IoSlice<'_>]) -> ::std::io::Result<usize> {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.write_vectored(bufs))
        }

        #[inline]
        fn write_all(&mut self, buf: &[u8]) -> ::std::io::Result<()> {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.write_all(buf))
        }

        #[inline]
        fn write_fmt(&mut self, fmt: ::std::fmt::Arguments<'_>) -> ::std::io::Result<()> {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.write_fmt(fmt))
        }
    };

    (@$dispatch:tt Iterator) => {
        #[inline]
        fn next(&mut self) -> Option<Self::Item> {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.next())
        }

        #[inline]
        fn size_hint(&self) -> (usize, Option<usize>) {
            $crate::impl_traits!(@$dispatch self, [&], inner => inner.size_hint())
        }

        #[inline]
        fn count(self) -> usize {
            $crate::impl_traits!(@$dispatch self, [], inner => inner.count())
        }

        #[inline]
        fn last(self) -> Option<Self::Item> {
            $crate::impl_traits!(@$dispatch self, [], inner => inner.last())
        }

        #[inline]
        fn nth(&mut self, n: usize) -> Option<Self::Item> {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.nth(n))
        }

        #[inline]
        fn for_each<F>(self, f: F)
        where
            Self: Sized,
            F: FnMut(Self::Item),
        {
            $crate::impl_traits!(@$dispatch self, [], inner => inner.for_each(f))
        }

        #[inline]
        fn collect<B>(self) -> B
        where
            B: FromIterator<Self::Item>,
            Self: Sized,
        {
            $crate::impl_traits!(@$dispatch self, [], inner => inner.collect())
        }

        #[inline]
        fn partition<B, F>(self, f: F) -> (B, B)
        where
            Self: Sized,
            B: Default + Extend<Self::Item>,
            F: FnMut(&Self::Item) -> bool,
        {
            $crate::impl_traits!(@$dispatch self, [], inner => inner.partition(f))
        }

        #[inline]
        fn try_fold<B, F, R>(&mut self, init: B, f: F) -> R
        where
            Self: Sized,
            F: FnMut(B, Self::Item) -> R,
            R: ::std::ops::Try<Output = B>,
        {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.try_fold(init, f))
        }

        #[inline]
        fn try_for_each<F, R>(&mut self, f: F) -> R
        where
            Self: Sized,
            F: FnMut(Self::Item) -> R,
            R: ::std::ops::Try<Output = ()>,
        {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.try_for_each(f))
        }

        #[inline]
        fn fold<B, F>(self, init: B, f: F) -> B
        where
            Self: Sized,
            F: FnMut(B, Self::Item) -> B,
        {
            $crate::impl_traits!(@$dispatch self, [], inner => inner.fold(init, f))
        }

        #[inline]
        fn reduce<F>(self, f: F) -> Option<Self::Item>
        where
            Self: Sized,
            F: FnMut(Self::Item, Self::Item) -> Self::Item,
        {
            $crate::impl_traits!(@$dispatch self, [], inner => inner.reduce(f))
        }

        #[inline]
        fn all<F>(&mut self, f: F) -> bool
        where
            Self: Sized,
            F: FnMut(Self::Item) -> bool,
        {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.all(f))
        }

        #[inline]
        fn any<F>(&mut self, f: F) -> bool
        where
            Self: Sized,
            F: FnMut(Self::Item) -> bool,
        {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.any(f))
        }

        #[inline]
        fn find<P>(&mut self, predicate: P) -> Option<Self::Item>
        where
            Self: Sized,
            P: FnMut(&Self::Item) -> bool,
        {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.find(predicate))
        }

        #[inline]
        fn find_map<B, F>(&mut self, f: F) -> Option<B>
        where
            Self: Sized,
            F: FnMut(Self::Item) -> Option<B>,
        {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.find_map(f))
        }

        #[inline]
        fn position<P>(&mut self, predicate: P) -> Option<usize>
        where
            Self: Sized,
            P: FnMut(Self::Item) -> bool,
        {
            $crate::impl_traits!(@$dispatch self, [&mut], inner => inner.position(predicate))
        }

        #[inline]
        fn max_by_key<B, F>(self, f: F) -> Option<Self::Item>
        where
            B: Ord,
            Self: Sized,
            F: FnMut(&Self::Item) -> B,
        {
            $crate::impl_traits!(@$dispatch self, [], inner => inner.max_by_key(f))
        }

        #[inline]
        fn max_by<F>(self, compare: F) -> Option<Self::Item>
        where
            Self: Sized,
            F: FnMut(&Self::Item, &Self::Item) -> ::std::cmp::Ordering,
        {
            $crate::impl_traits!(@$dispatch self, [], inner => inner.max_by(compare))
        }

        #[inline]
        fn min_by_key<B, F>(self, f: F) -> Option<Self::Item>
        where
            B: Ord,
            Self: Sized,
            F: FnMut(&Self::Item) -> B,
        {
            $crate::impl_traits!(@$dispatch self, [], inner => inner.min_by_key(f))
        }

        #[inline]
        fn min_by<F>(self, compare: F) -> Option<Self::Item>
        where
            Self: Sized,
            F: FnMut(&Self::Item, &Self::Item) -> ::std::cmp::Ordering,
        {
            $crate::impl_traits!(@$dispatch self, [], inner => inner.min_by(compare))
        }

        #[inline]
        fn sum<S>(self) -> S
        where
            Self: Sized,
            S: ::std::iter::Sum<Self::Item>,
        {
            $crate::impl_traits!(@$dispatch self, [], inner => inner.sum())
        }

        #[inline]
        fn product<S>(self) -> S
        where
            Self: Sized,
            S: ::std::iter::Product<Self::Item>,
        {
            $crate::impl_traits!(@$dispatch self, [], inner => inner.product())
        }

        #[inline]
        fn is_sorted_by<F>(self, compare: F) -> bool
        where
            Self: Sized,
            F: FnMut(&Self::Item, &Self::Item) -> bool,
        {
            $crate::impl_traits!(@$dispatch self, [], inner => inner.is_sorted_by(compare))
        }

        #[inline]
        fn is_sorted_by_key<F, K>(self, f: F) -> bool
        where
            Self: Sized,
            F: FnMut(Self::Item) -> K,
            K: PartialOrd,
        {
            $crate::impl_traits!(@$dispatch self, [], inner => inner.is_sorted_by_key(f))
        }
    };

    (@match_stmt $value:expr, [$($reference:tt)*], $pattern:pat => $result:expr) => {
        match_macro!($value, $pattern => $result)
    };

    (@wrapper_type $value:expr, [$($reference:tt)*], $pattern:pat => $result:expr) => {
        {
            let $pattern = $($reference)* $value.0;
            $result
        }
    };
}
