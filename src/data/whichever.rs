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
///     /// An enum representing the allowable input types
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
///     /// An enum representing iterators yielding the same value several times
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
/// Any trait being implemented must be imported (no qualified paths are
/// possible with this macro).
///
/// When implementing [`Iterator`], the generics `Ff`, `Bb`, `Rr`, `Pp`, `Ss`,
/// and `Kk` are reserved. These were chosen because they are unlikely to
/// conflict with generics on the type. Also, the `try_trait_v2` feature must be
/// enabled.
///
/// For [`Iterator`], it may also be desirable to have enum variants whose
/// `Item` associated types differ, but use some map operation to unify them.
/// This is possible with the following syntax:
///
/// ```
/// #![feature(try_trait_v2)]
/// use zoe::define_whichever;
///
/// define_whichever!{
///     pub enum SignedOrUnsigned {
///         Signed(std::iter::RepeatN<i32>),
///         Unsigned(std::iter::RepeatN<u32>),
///     }
///
///     #[map(i64::from)]
///     impl Iterator for SignedOrUnsigned {
///         type Item = i64;
///     }
/// }
/// ```
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

        $crate::impl_traits! {
            @[match_enum $struct_name { $($variant),+ }]
            $($impl_stuff)*
        }
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
/// Any trait being implemented must be imported (no qualified paths are
/// possible with this macro).
///
/// When implementing [`Iterator`], the generics `Ff`, `Bb`, `Rr`, `Pp`, `Ss`,
/// and `Kk` are reserved. These were chosen because they are unlikely to
/// conflict with generics on the type. Also, the `try_trait_v2` feature must be
/// enabled.
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
        $crate::impl_traits! {
            @[wrapper_type]
            impl $($tt)*
        }
    };

    (@$dispatch:tt
        $(
            $(#[map($map:expr)])?
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

                $crate::impl_traits! {
                    @methods [$dispatch $($map)?] $trait
                }
            }
        )*
    };

    (@methods [$dispatch:tt] Read) => {
        #[inline]
        fn read(&mut self, buf: &mut [u8]) -> ::std::io::Result<usize> {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner.read(buf))
        }

        #[inline]
        fn read_vectored(&mut self, bufs: &mut [::std::io::IoSliceMut<'_>]) -> ::std::io::Result<usize> {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner.read_vectored(bufs))
        }

        #[inline]
        fn read_to_end(&mut self, buf: &mut ::std::vec::Vec<u8>) -> ::std::io::Result<usize> {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner.read_to_end(buf))
        }

        #[inline]
        fn read_to_string(&mut self, buf: &mut ::std::string::String) -> ::std::io::Result<usize> {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner.read_to_string(buf))
        }

        #[inline]
        fn read_exact(&mut self, buf: &mut [u8]) -> ::std::io::Result<()> {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner.read_exact(buf))
        }
    };

    (@methods [$dispatch:tt] Write) => {
        #[inline]
        fn write(&mut self, buf: &[u8]) -> ::std::io::Result<usize> {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner.write(buf))
        }

        #[inline]
        fn flush(&mut self) -> ::std::io::Result<()> {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner.flush())
        }

        #[inline]
        fn write_vectored(&mut self, bufs: &[::std::io::IoSlice<'_>]) -> ::std::io::Result<usize> {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner.write_vectored(bufs))
        }

        #[inline]
        fn write_all(&mut self, buf: &[u8]) -> ::std::io::Result<()> {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner.write_all(buf))
        }

        #[inline]
        fn write_fmt(&mut self, fmt: ::std::fmt::Arguments<'_>) -> ::std::io::Result<()> {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner.write_fmt(fmt))
        }
    };

    (@methods [$dispatch:tt $($map:expr)?] Iterator) => {
        #[inline]
        fn next(&mut self) -> Option<Self::Item> {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner.next()$(.map($map))?)
        }

        #[inline]
        fn size_hint(&self) -> (usize, Option<usize>) {
            $crate::impl_traits!(@delegate $dispatch self, [&], inner => inner.size_hint())
        }

        #[inline]
        fn count(self) -> usize {
            $crate::impl_traits!(@delegate $dispatch self, [], inner => inner.count())
        }

        #[inline]
        fn last(self) -> Option<Self::Item> {
            $crate::impl_traits!(@delegate $dispatch self, [], inner => inner.last()$(.map($map))?)
        }

        #[inline]
        fn nth(&mut self, n: usize) -> Option<Self::Item> {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner.nth(n)$(.map($map))?)
        }

        #[inline]
        fn for_each<Ff>(self, f: Ff)
        where
            Self: Sized,
            Ff: FnMut(Self::Item),
        {
            $crate::impl_traits!(@delegate $dispatch self, [], inner => inner$(.map($map))?.for_each(f))
        }

        #[inline]
        fn collect<Bb>(self) -> Bb
        where
            Bb: FromIterator<Self::Item>,
            Self: Sized,
        {
            $crate::impl_traits!(@delegate $dispatch self, [], inner => inner$(.map($map))?.collect())
        }

        #[inline]
        fn partition<Bb, Ff>(self, f: Ff) -> (Bb, Bb)
        where
            Self: Sized,
            Bb: Default + Extend<Self::Item>,
            Ff: FnMut(&Self::Item) -> bool,
        {
            $crate::impl_traits!(@delegate $dispatch self, [], inner => inner$(.map($map))?.partition(f))
        }

        #[inline]
        fn try_fold<Bb, Ff, Rr>(&mut self, init: Bb, f: Ff) -> Rr
        where
            Self: Sized,
            Ff: FnMut(Bb, Self::Item) -> Rr,
            Rr: ::std::ops::Try<Output = Bb>,
        {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner$(.map($map))?.try_fold(init, f))
        }

        #[inline]
        fn try_for_each<Ff, Rr>(&mut self, f: Ff) -> Rr
        where
            Self: Sized,
            Ff: FnMut(Self::Item) -> Rr,
            Rr: ::std::ops::Try<Output = ()>,
        {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner$(.map($map))?.try_for_each(f))
        }

        #[inline]
        fn fold<Bb, Ff>(self, init: Bb, f: Ff) -> Bb
        where
            Self: Sized,
            Ff: FnMut(Bb, Self::Item) -> Bb,
        {
            $crate::impl_traits!(@delegate $dispatch self, [], inner => inner$(.map($map))?.fold(init, f))
        }

        #[inline]
        fn reduce<Ff>(self, f: Ff) -> Option<Self::Item>
        where
            Self: Sized,
            Ff: FnMut(Self::Item, Self::Item) -> Self::Item,
        {
            $crate::impl_traits!(@delegate $dispatch self, [], inner => inner$(.map($map))?.reduce(f))
        }

        #[inline]
        fn all<Ff>(&mut self, f: Ff) -> bool
        where
            Self: Sized,
            Ff: FnMut(Self::Item) -> bool,
        {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner$(.map($map))?.all(f))
        }

        #[inline]
        fn any<Ff>(&mut self, f: Ff) -> bool
        where
            Self: Sized,
            Ff: FnMut(Self::Item) -> bool,
        {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner$(.map($map))?.any(f))
        }

        #[inline]
        fn find<Pp>(&mut self, predicate: Pp) -> Option<Self::Item>
        where
            Self: Sized,
            Pp: FnMut(&Self::Item) -> bool,
        {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner$(.map($map))?.find(predicate))
        }

        #[inline]
        fn find_map<Bb, Ff>(&mut self, f: Ff) -> Option<Bb>
        where
            Self: Sized,
            Ff: FnMut(Self::Item) -> Option<Bb>,
        {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner$(.map($map))?.find_map(f))
        }

        #[inline]
        fn position<Pp>(&mut self, predicate: Pp) -> Option<usize>
        where
            Self: Sized,
            Pp: FnMut(Self::Item) -> bool,
        {
            $crate::impl_traits!(@delegate $dispatch self, [&mut], inner => inner$(.map($map))?.position(predicate))
        }

        #[inline]
        fn max_by_key<Bb, Ff>(self, f: Ff) -> Option<Self::Item>
        where
            Bb: Ord,
            Self: Sized,
            Ff: FnMut(&Self::Item) -> Bb,
        {
            $crate::impl_traits!(@delegate $dispatch self, [], inner => inner$(.map($map))?.max_by_key(f))
        }

        #[inline]
        fn max_by<Ff>(self, compare: Ff) -> Option<Self::Item>
        where
            Self: Sized,
            Ff: FnMut(&Self::Item, &Self::Item) -> ::std::cmp::Ordering,
        {
            $crate::impl_traits!(@delegate $dispatch self, [], inner => inner$(.map($map))?.max_by(compare))
        }

        #[inline]
        fn min_by_key<Bb, Ff>(self, f: Ff) -> Option<Self::Item>
        where
            Bb: Ord,
            Self: Sized,
            Ff: FnMut(&Self::Item) -> Bb,
        {
            $crate::impl_traits!(@delegate $dispatch self, [], inner => inner$(.map($map))?.min_by_key(f))
        }

        #[inline]
        fn min_by<Ff>(self, compare: Ff) -> Option<Self::Item>
        where
            Self: Sized,
            Ff: FnMut(&Self::Item, &Self::Item) -> ::std::cmp::Ordering,
        {
            $crate::impl_traits!(@delegate $dispatch self, [], inner => inner$(.map($map))?.min_by(compare))
        }

        #[inline]
        fn sum<Ss>(self) -> Ss
        where
            Self: Sized,
            Ss: ::std::iter::Sum<Self::Item>,
        {
            $crate::impl_traits!(@delegate $dispatch self, [], inner => inner$(.map($map))?.sum())
        }

        #[inline]
        fn product<Ss>(self) -> Ss
        where
            Self: Sized,
            Ss: ::std::iter::Product<Self::Item>,
        {
            $crate::impl_traits!(@delegate $dispatch self, [], inner => inner$(.map($map))?.product())
        }

        #[inline]
        fn is_sorted_by<Ff>(self, compare: Ff) -> bool
        where
            Self: Sized,
            Ff: FnMut(&Self::Item, &Self::Item) -> bool,
        {
            $crate::impl_traits!(@delegate $dispatch self, [], inner => inner$(.map($map))?.is_sorted_by(compare))
        }

        #[inline]
        fn is_sorted_by_key<Ff, Kk>(self, f: Ff) -> bool
        where
            Self: Sized,
            Ff: FnMut(Self::Item) -> Kk,
            Kk: PartialOrd,
        {
            $crate::impl_traits!(@delegate $dispatch self, [], inner => inner$(.map($map))?.is_sorted_by_key(f))
        }
    };

    (@delegate [match_enum $enum_name:ident { $($variant:ident),+ }] $value:expr, [$($reference:tt)*], $pattern:pat => $result:expr) => {
        match $value {
            $(
                $enum_name::$variant($pattern) => $result,
            )+
        }
    };

    (@delegate [wrapper_type] $value:expr, [$($reference:tt)*], $pattern:pat => $result:expr) => {
        {
            let $pattern = $($reference)* $value.0;
            $result
        }
    };
}
