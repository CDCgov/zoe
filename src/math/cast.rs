// TODO: If non-lifetime binders with where clauses become available, then we
// can use implication predicates (item 3 of
// https://github.com/rust-lang/types-team/issues/81) to write a trait bound
// requiring AnyInt to be able to be cast to AnyInt. Once this becomes
// available, Cast and CastFrom can be combined and generalized. Or, if
// elaboration of where clauses becomes available, we can also simplify the
// machinery here (https://github.com/rust-lang/rust/issues/20671).

mod private {
    /// A sealed helper trait for specifying pairs of primitives where T can be
    /// cast to S. If such a cast is possible, then [`CastPair`] is implemented
    /// for `()`.
    pub trait CastPair<T, S> {
        /// Casts `n` of type `T` to type `S`.
        fn cast(n: T) -> S;
    }

    /// Implements [`CastPair`] for primitive types.
    macro_rules! impl_cast_pairs {
        ($($from:ty => $($to:ty),*;)*) => {
            $(
                $(
                    impl CastPair<$from, $to> for () {
                        #[allow(clippy::cast_lossless)]
                        #[allow(clippy::cast_sign_loss)]
                        #[allow(clippy::cast_possible_wrap)]
                        #[allow(clippy::cast_precision_loss)]
                        #[allow(clippy::cast_possible_truncation)]
                        fn cast(n: $from) -> $to {
                            n as $to
                        }
                    }
                )*
            )*
        };
    }

    // No type can be cast to a bool other than itself, and bools cannot be cast
    // to floats
    impl_cast_pairs! {
        u8 => u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64;
        u16 => u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64;
        u32 => u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64;
        u64 => u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64;
        u128 => u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64;
        usize => u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64;
        i8 => u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64;
        i16 => u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64;
        i32 => u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64;
        i64 => u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64;
        i128 => u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64;
        isize => u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64;
        f32 => u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64;
        f64 => u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64;
        bool => bool, u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize;
    }

    /// A helper trait for [`CastAsFixed`]. In order to the trait bounds to be
    /// elaborated (meaning that where clauses for [`CastPair`] do not need to
    /// be added any time we cast types), the trait bounds must be bounds on
    /// `Self`. This trait changes the first generic in [`CastPair`] to `Self`,
    /// and uses associated type `P` along with the equality constraint in
    /// [`CastAsFixed`] to ensure [`CastPair`] is implemented for `()`.
    pub trait CastAsHelper<S>: Sized {
        type P: CastPair<Self, S>;

        /// Casts `self` to the type `S` specified in the trait.
        fn cast_as_helper(self) -> S {
            Self::P::cast(self)
        }
    }

    impl<T, S> CastAsHelper<S> for T
    where
        (): CastPair<T, S>,
    {
        type P = ();
    }

    /// A trait allowing casting of `Self` into `S`. Adding this trait as a
    /// bound on another trait ensures that [`CastAs`] can be used on types
    /// implementing that trait.
    ///
    /// [`CastAs`]: super::CastAs
    pub trait CastAsFixed<S>: CastAsHelper<S, P = ()> {
        /// Casts `self` to the type `S` specified in the trait.
        fn cast_as_fixed(self) -> S {
            <Self as CastAsHelper<S>>::cast_as_helper(self)
        }
    }

    impl<T: CastAsHelper<S, P = ()>, S> CastAsFixed<S> for T {}

    /// A helper trait for [`CastFromFixed`]. In order to the trait bounds to be
    /// elaborated (meaning that where clauses for [`CastPair`] do not need to
    /// be added any time we cast types), the trait bounds must be bounds on
    /// `Self`. This trait changes the first generic in [`CastPair`] to `Self`,
    /// and uses associated type `P` along with the equality constraint in
    /// [`CastFromFixed`] to ensure [`CastPair`] is implemented for `()`.
    pub trait CastFromHelper<T>: Sized {
        type P: CastPair<T, Self>;

        /// Casts `val` from the type `T` specified in the trait to `Self`.
        fn cast_from_helper(val: T) -> Self {
            Self::P::cast(val)
        }
    }

    impl<T, S> CastFromHelper<T> for S
    where
        (): CastPair<T, S>,
    {
        type P = ();
    }

    /// A trait allowing casting of `T` into `Self`. Adding this trait as a
    /// bound on another trait ensures that [`CastFrom`] can be used on types
    /// implementing that trait.
    ///
    /// [`CastFrom`]: super::CastFrom
    pub trait CastFromFixed<T>: CastFromHelper<T, P = ()> {
        /// Casts `val` from the type `T` specified in the trait to `Self`.
        fn cast_from_fixed(val: T) -> Self {
            <Self as CastFromHelper<T>>::cast_from_helper(val)
        }
    }

    impl<T, S: CastFromHelper<T, P = ()>> CastFromFixed<T> for S {}
}

/// A trait providing the ability to cast from a potentially generic type (such
/// as [`AnyInt`]) to a known primitive type.
///
/// [`AnyInt`]: crate::math::AnyInt
pub trait CastAs {
    /// Performs a primitive cast from `self` (a potentially generic type such
    /// as [`AnyInt`]) to a known primitive type.
    ///
    /// ## Examples
    ///
    /// ```
    /// # use zoe::math::AnyInt;
    /// fn as_f64<T: AnyInt>(val: T) -> f64 {
    ///     val.cast_as::<f64>()
    /// }
    /// ```
    ///
    /// [`AnyInt`]: crate::math::AnyInt
    #[inline]
    #[must_use]
    fn cast_as<S>(self) -> S
    where
        Self: private::CastAsFixed<S>, {
        <Self as private::CastAsFixed<S>>::cast_as_fixed(self)
    }
}

impl<T: Primitive> CastAs for T {}

/// A trait providing the ability to cast from a known primitive type to a
/// potentially generic type (such as [`AnyInt`]).
///
/// [`AnyInt`]: crate::math::AnyInt
pub trait CastFrom: Sized {
    /// Performs a primitive cast from `val` (a known primitive type) to a
    /// potentially generic type such as [`AnyInt`].
    ///
    /// ## Examples
    ///
    /// ```
    /// # use zoe::math::AnyInt;
    /// fn from_f64<T: AnyInt>(val: f64) -> T {
    ///     T::cast_from(val)
    /// }
    /// ```
    ///
    /// [`AnyInt`]: crate::math::AnyInt
    fn cast_from<T>(val: T) -> Self
    where
        Self: private::CastFromFixed<T>, {
        <Self as private::CastFromFixed<T>>::cast_from_fixed(val)
    }
}

impl<T: Primitive> CastFrom for T {}

/// A trait which ensures `Self` can be cast to any known (non-generic)
/// numeric primitive.
pub trait CastAsNumeric:
    CastAs
    + private::CastAsFixed<u8>
    + private::CastAsFixed<u16>
    + private::CastAsFixed<u32>
    + private::CastAsFixed<u64>
    + private::CastAsFixed<u128>
    + private::CastAsFixed<usize>
    + private::CastAsFixed<i8>
    + private::CastAsFixed<i32>
    + private::CastAsFixed<i64>
    + private::CastAsFixed<i128>
    + private::CastAsFixed<isize>
    + private::CastAsFixed<f32>
    + private::CastAsFixed<f64> {
}

impl<T> CastAsNumeric for T where
    T: CastAs
        + private::CastAsFixed<u8>
        + private::CastAsFixed<u16>
        + private::CastAsFixed<u32>
        + private::CastAsFixed<u64>
        + private::CastAsFixed<u128>
        + private::CastAsFixed<usize>
        + private::CastAsFixed<i8>
        + private::CastAsFixed<i32>
        + private::CastAsFixed<i64>
        + private::CastAsFixed<i128>
        + private::CastAsFixed<isize>
        + private::CastAsFixed<f32>
        + private::CastAsFixed<f64>
{
}

/// A trait which ensures any known (non-generic) numeric primitive can be cast
/// to `Self`.
pub trait CastFromNumeric:
    CastFrom
    + private::CastFromFixed<u8>
    + private::CastFromFixed<u16>
    + private::CastFromFixed<u32>
    + private::CastFromFixed<u64>
    + private::CastFromFixed<u128>
    + private::CastFromFixed<usize>
    + private::CastFromFixed<i8>
    + private::CastFromFixed<i32>
    + private::CastFromFixed<i64>
    + private::CastFromFixed<i128>
    + private::CastFromFixed<isize>
    + private::CastFromFixed<f32>
    + private::CastFromFixed<f64> {
}

impl<T> CastFromNumeric for T where
    T: CastFrom
        + private::CastFromFixed<u8>
        + private::CastFromFixed<u16>
        + private::CastFromFixed<u32>
        + private::CastFromFixed<u64>
        + private::CastFromFixed<u128>
        + private::CastFromFixed<usize>
        + private::CastFromFixed<i8>
        + private::CastFromFixed<i32>
        + private::CastFromFixed<i64>
        + private::CastFromFixed<i128>
        + private::CastFromFixed<isize>
        + private::CastFromFixed<f32>
        + private::CastFromFixed<f64>
{
}

/// A trait which ensures that a boolean can be cast to `Self`.
pub trait CastFromBool: private::CastFromFixed<bool> {}

impl<T: private::CastFromFixed<bool>> CastFromBool for T {}

/// A marker trait for primitive numeric types or booleans.
pub trait Primitive: Sized {}

impl Primitive for u8 {}
impl Primitive for u16 {}
impl Primitive for u32 {}
impl Primitive for u64 {}
impl Primitive for u128 {}
impl Primitive for usize {}
impl Primitive for i8 {}
impl Primitive for i16 {}
impl Primitive for i32 {}
impl Primitive for i64 {}
impl Primitive for i128 {}
impl Primitive for isize {}
impl Primitive for f32 {}
impl Primitive for f64 {}
impl Primitive for bool {}
