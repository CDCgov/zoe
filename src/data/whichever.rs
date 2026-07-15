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
/// possible with this macro). Const generics are not supported on the type.
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
/// [`Read`]: std::io::Read
/// [`Write`]: std::io::Write
#[macro_export]
macro_rules! define_whichever {
    // The public entry-point for the macro
    (
        $(#[$meta:meta])*
        $vis:vis enum $struct_name:ident $($tt:tt)*
    ) => {
        $crate::define_whichever! {
            @parse_header
            [$(#[$meta])*] [$vis] [$struct_name] []
            $($tt)*
        }
    };

    // parse_header is a tt-muncher which attempts to separate generics and
    // where clauses on the enum from the variants and impl statements

    // The variants have been found, so emit the enum and delegate to
    // impl_traits
    (@parse_header
        [$(#[$meta:meta])*] [$vis:vis] [$struct_name:ident] [$($header:tt)*]
        {
            $(
                $(#[$variant_meta:meta])*
                $variant:ident($ty:ty)
            ),+
            $(,)?
        }
        $($impl_stuff:tt)*
    ) => {
        $(#[$meta])*
        $vis enum $struct_name $($header)* {
            $(
                $(#[$variant_meta])*
                $variant($ty),
            )+
        }

        $crate::impl_traits! {
            @parse_impls [match_enum $struct_name { $($variant),+ }]
            $($impl_stuff)*
        }
    };

    // Other case: continue munching
    (@parse_header
        [$($meta:tt)*] [$vis:vis] [$struct_name:ident] [$($header:tt)*]
        $tt:tt $($rest:tt)*
    ) => {
        $crate::define_whichever! {
            @parse_header
            [$($meta)*] [$vis] [$struct_name] [$($header)* $tt]
            $($rest)*
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
/// possible with this macro). Const generics are not supported on the type.
///
/// When implementing [`Iterator`], the generics `Ff`, `Bb`, `Rr`, `Pp`, `Ss`,
/// and `Kk` are reserved. These were chosen because they are unlikely to
/// conflict with generics on the type. Also, the `try_trait_v2` feature must be
/// enabled.
///
/// [`Read`]: std::io::Read
/// [`Write`]: std::io::Write
#[macro_export]
macro_rules! impl_traits {
    // The public entry-point for the macro, which assumes that the traits are
    // being implemented on a wrapper type
    (impl $($tt:tt)*) => {
        $crate::impl_traits! {
            @parse_impls [wrapper_type]
            impl $($tt)*
        }
    };

    // parse_impls is the first step in parsing a single impl statement. When
    // that impl statement finishes parsing, it returns back to the
    // `parse_impls` branch to restart again on the next one.

    // When no more impls remain, we are done
    (@parse_impls $dispatch:tt) => {};

    // Parse the header of the impl statement
    (@parse_impls $dispatch:tt
        $(#[map($map:expr)])?
        impl $($rest:tt)*
    ) => {
        $crate::impl_traits! {
            @parse_impl_header
            $dispatch [$($map)?] []
            $($rest)*
        }
    };

    // While parse_impl_header is specified, the macro attempts to parse the
    // generics on impl statement, before the trait name. Tokens are munched
    // until `trait for type where ...` is found or `trait for type { ... }` is
    // found.

    // Case where `trait for type where ...` is found
    (@parse_impl_header
        $dispatch:tt [$($map:expr)?] [$($header:tt)*]
        $trait:ident for $self_ty:ty where $($rest:tt)*
    ) => {
        $crate::impl_traits! {
            @parse_where
            $dispatch [$($map)?] [$($header)*] [$trait] [$self_ty] []
            $($rest)*
        }
    };

    // Case where `trait for type { ... }` is found (no where clause)
    (@parse_impl_header
        $dispatch:tt [$($map:expr)?] [$($header:tt)*]
        $trait:ident for $self_ty:ty {
            $($body:tt)*
        }
        $($rest:tt)*
    ) => {
        $crate::impl_traits! {
            @emit_impl
            $dispatch [$($map)?] [$($header)*] [$trait] [$self_ty] []
            { $($body)* }
        }
        $crate::impl_traits! {
            @parse_impls $dispatch $($rest)*
        }
    };

    // Other case: continue munching
    (@parse_impl_header
        $dispatch:tt [$($map:expr)?] [$($header:tt)*]
        $tt:tt $($rest:tt)*
    ) => {
        $crate::impl_traits! {
            @parse_impl_header
            $dispatch [$($map)?] [$($header)* $tt]
            $($rest)*
        }
    };

    // While parse_where is specified, the macro attempts to parse the where
    // clause. Tokens are munched until `{ ... }` is found, either ending the
    // arguments or immediately followed by another `impl` statement,
    // potentially annotated with `map`.

    // Case where `{ ... }` is next
    (@parse_where
        $dispatch:tt [$($map:expr)?] [$($header:tt)*]
        [$trait:ident] [$self_ty:ty] [$($where_clause:tt)*]
        {
            $($body:tt)*
        }
        $($(#[map($next_map:expr)])? impl $($rest:tt)*)?
    ) => {
        $crate::impl_traits! {
            @finish_where
            $dispatch [$($map)?] [$($header)*] [$trait] [$self_ty]
            [$($where_clause)*] { $($body)* }
            [$($(#[map($next_map)])? impl $($rest)*)?]
        }
    };

    // Other case: continue munching
    (@parse_where
        $dispatch:tt [$($map:expr)?] [$($header:tt)*]
        [$trait:ident] [$self_ty:ty] [$($where_clause:tt)*]
        $tt:tt $($rest:tt)*
    ) => {
        $crate::impl_traits! {
            @parse_where
            $dispatch [$($map)?] [$($header)*]
            [$trait] [$self_ty] [$($where_clause)* $tt]
            $($rest)*
        }
    };

    // Finish a where clause by emitting the impl and then continuing to parse
    // the next one
    (@finish_where
        $dispatch:tt [$($map:expr)?] [$($header:tt)*]
        [$trait:ident] [$self_ty:ty] [$($where_clause:tt)*]
        {
            $($body:tt)*
        }
        [$($rest:tt)*]
    ) => {
        $crate::impl_traits! {
            @emit_impl
            $dispatch [$($map)?] [$($header)*] [$trait] [$self_ty]
            [where $($where_clause)*]
            { $($body)* }
        }
        $crate::impl_traits! {
            @parse_impls $dispatch $($rest)*
        }
    };

    // Emit a single impl statement, given that the tokens have been grouped
    // together by the tt-muncher
    (@emit_impl
        $dispatch:tt [$($map:expr)?] [$($header:tt)*]
        [$trait:ident] [$self_ty:ty] [$($where_clause:tt)*]
        {
            $($body:tt)*
        }
    ) => {
        impl $($header)* $trait for $self_ty $($where_clause)* {
            $($body)*

            $crate::impl_traits! {
                @methods [$dispatch $($map)?] $trait
            }
        }
    };

    // Generates the methods for the Read trait
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

    // Generates the methods for the Write trait
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

    // Generates the methods for the Iterator trait
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

    // Provides the body for a single method implementation, given that a match
    // statement is being used (this path is taken by define_whichever)
    (@delegate [match_enum $enum_name:ident { $($variant:ident),+ }] $value:expr, [$($reference:tt)*], $pattern:pat => $result:expr) => {
        match $value {
            $(
                $enum_name::$variant($pattern) => $result,
            )+
        }
    };

    // Provides the body for a single method implementation, given that a
    // wrapper type is used (this path is taken by the public entrypoint of
    // impl_traits)
    (@delegate [wrapper_type] $value:expr, [$($reference:tt)*], $pattern:pat => $result:expr) => {
        {
            let $pattern = $($reference)* $value.0;
            $result
        }
    };
}

#[cfg(test)]
mod tests {
    #[test]
    fn wrapper_supports_lifetimes_and_type_parameters() {
        struct MyIter<'a, T>(std::slice::Iter<'a, T>);

        impl_traits! {
            impl<'a, T> Iterator for MyIter<'a, T> {
                type Item = &'a T;
            }
        }

        let values = [1, 2, 3];
        let mut iter = MyIter(values.iter());

        assert_eq!(iter.next(), Some(&1));
        assert_eq!(iter.collect::<Vec<_>>(), vec![&2, &3]);
    }

    #[test]
    fn enum_supports_lifetimes_and_type_parameters() {
        define_whichever! {
            enum MyIter<'a, T> {
                Slice(std::slice::Iter<'a, T>),
                Once(std::iter::Once<&'a T>),
            }

            impl<'a, T> Iterator for MyIter<'a, T> {
                type Item = &'a T;
            }
        }

        let values = [1, 2, 3];
        let mut slice = MyIter::Slice(values.iter());
        let mut once = MyIter::Once(std::iter::once(&values[0]));

        assert_eq!(slice.next(), Some(&1));
        assert_eq!(once.next(), Some(&1));
    }

    #[test]
    fn wrapper_supports_const_parameters() {
        struct MyIter<const N: usize>(std::array::IntoIter<u8, N>);

        impl_traits! {
            impl<const N: usize> Iterator for MyIter<N> {
                type Item = u8;
            }
        }

        let iter = MyIter([1, 2, 3].into_iter());

        assert_eq!(iter.collect::<Vec<_>>(), vec![1, 2, 3]);
    }

    #[test]
    fn enum_supports_const_parameters() {
        define_whichever! {
            enum MyIter<const N: usize> {
                Array(std::array::IntoIter<u8, N>),
            }

            impl<const N: usize> Iterator for MyIter<N> {
                type Item = u8;
            }
        }

        let iter = MyIter::Array([1, 2, 3].into_iter());

        assert_eq!(iter.collect::<Vec<_>>(), vec![1, 2, 3]);
    }

    #[test]
    fn wrapper_supports_mixed_parameters() {
        pub struct MyIter<'a, T, const N: usize>(std::slice::Iter<'a, [T; N]>);

        impl_traits! {
            impl<'a, T, const N: usize> Iterator for MyIter<'a, T, N> {
                type Item = &'a [T; N];
            }
        }

        let values = [[1, 2], [3, 4]];
        let mut iter = MyIter(values.iter());

        assert_eq!(iter.next(), Some(&[1, 2]));
        assert_eq!(iter.next(), Some(&[3, 4]));
    }

    #[test]
    fn enum_supports_mixed_parameters() {
        define_whichever! {
            enum MyIter<'a, T, const N: usize> {
                Iter(std::slice::Iter<'a, [T; N]>),
            }

            impl<'a, T, const N: usize> Iterator for MyIter<'a, T, N> {
                type Item = &'a [T; N];
            }
        }

        let values = [[1, 2], [3, 4]];
        let mut iter = MyIter::Iter(values.iter());

        assert_eq!(iter.next(), Some(&[1, 2]));
        assert_eq!(iter.next(), Some(&[3, 4]));
    }

    #[test]
    fn wrapper_supports_multiple_where_bounds() {
        pub struct MyIter<T>(std::vec::IntoIter<T>);

        impl_traits! {
            impl<T: Clone> Iterator for MyIter<T> where T: Sync + Send {
                type Item = T;
            }
        }

        let iter = MyIter(vec![1, 2, 3].into_iter());

        assert_eq!(iter.collect::<Vec<_>>(), vec![1, 2, 3]);
    }

    #[test]
    fn enum_supports_multiple_where_bounds() {
        define_whichever! {
            enum MyIter<T> {
                Iter(std::vec::IntoIter<T>)
            }

            impl<T: Clone> Iterator for MyIter<T> where T: Sync + Send {
                type Item = T;
            }
        }

        let iter = MyIter::Iter(vec![1, 2, 3].into_iter());

        assert_eq!(iter.collect::<Vec<_>>(), vec![1, 2, 3]);
    }
}
