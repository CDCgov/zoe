#[cfg(feature = "fuzzing")]
pub mod float_compare {
    use crate::{
        alignment::phmm::{CorePhmm, EmissionParams, GlobalPhmm, LayerParams, TransitionParams},
        math::{NearlyEqual, NearlyEqualMethod},
    };

    impl<T: NearlyEqual<T> + Copy> NearlyEqual<T> for TransitionParams<T> {
        fn nearly_equal<M: NearlyEqualMethod<T>>(&self, b: &Self, strategy: &M) -> (bool, Option<(T, T)>) {
            self.0.nearly_equal(&b.0, strategy)
        }
    }

    impl<T: NearlyEqual<T> + Copy, const S: usize> NearlyEqual<T> for EmissionParams<T, S> {
        fn nearly_equal<M: NearlyEqualMethod<T>>(&self, b: &Self, strategy: &M) -> (bool, Option<(T, T)>) {
            self.0.nearly_equal(&b.0, strategy)
        }
    }

    impl<T: NearlyEqual<T> + Copy, const S: usize> NearlyEqual<T> for LayerParams<T, S> {
        fn nearly_equal<M: NearlyEqualMethod<T>>(&self, b: &Self, strategy: &M) -> (bool, Option<(T, T)>) {
            let (eq, vals) = self.transition.nearly_equal(&b.transition, strategy);
            if !eq {
                return (false, vals);
            }
            let (eq, vals) = self.emission_match.nearly_equal(&b.emission_match, strategy);
            if !eq {
                return (false, vals);
            }
            let (eq, vals) = self.emission_insert.nearly_equal(&b.emission_insert, strategy);
            if eq { (true, None) } else { (false, vals) }
        }
    }

    impl<T: NearlyEqual<T> + Copy, const S: usize> NearlyEqual<T> for CorePhmm<T, S> {
        fn nearly_equal<M: NearlyEqualMethod<T>>(&self, b: &Self, strategy: &M) -> (bool, Option<(T, T)>) {
            self.0.nearly_equal(&b.0, strategy)
        }
    }

    impl<T: NearlyEqual<T> + Copy, const S: usize> NearlyEqual<T> for GlobalPhmm<T, S> {
        fn nearly_equal<M: NearlyEqualMethod<T>>(&self, b: &Self, strategy: &M) -> (bool, Option<(T, T)>) {
            if self.mapping == b.mapping {
                self.core.nearly_equal(&b.core, strategy)
            } else {
                (false, None)
            }
        }
    }
}
