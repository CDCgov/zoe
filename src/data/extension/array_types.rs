// TODO: Make into a trait when functionality becomes available.
// See: https://github.com/rust-lang/rust-project-goals/issues/106

use crate::math::AnyInt;

/// Converts all bytes in an array to ASCII uppercase using
/// [`u8::to_ascii_uppercase`].
pub(crate) const fn make_uppercase<const N: usize>(a: &[u8; N]) -> [u8; N] {
    let mut b = [0; N];
    let mut i = 0;
    while i < a.len() {
        b[i] = a[i].to_ascii_uppercase();
        i += 1;
    }
    b
}

/// Identifies the first index of `needle` in `a`.
pub(crate) const fn position<const N: usize>(a: &[u8; N], needle: u8) -> Option<usize> {
    let mut i = 0;
    while i < a.len() {
        if needle == a[i] {
            return Some(i);
        }
        i += 1;
    }
    None
}

/// Returns whether all elements in `a` are unique.
pub(crate) const fn is_unique<const N: usize>(a: &[u8; N]) -> bool {
    let mut i = 0;
    while i < a.len() {
        let mut j = i + 1;
        while j < a.len() {
            if a[i] == a[j] {
                return false;
            }
            j += 1;
        }
        i += 1;
    }
    true
}

/// Checks whether all the values in `subset` also appear in `arr`.
pub(crate) const fn is_subset<const N: usize, const L: usize>(arr: &[u8; N], subset: &[u8; L]) -> bool {
    let mut i = 0;
    while i < subset.len() {
        let mut j = 0;
        let mut match_found = false;
        while j < arr.len() {
            if arr[j] == subset[i] {
                match_found = true;
                break;
            }
            j += 1;
        }
        if !match_found {
            return false;
        }
        i += 1;
    }
    true
}

// TODO: Update to std if it becomes available
// See: https://github.com/rust-lang/rust/issues/78504
#[allow(dead_code)]
pub(crate) const fn arr_max<const N: usize>(a: &[u8; N]) -> Option<u8> {
    let mut out = None;
    let mut i = 0;
    while i < a.len() {
        match out {
            Some(max) => {
                if a[i] > max {
                    out = Some(a[i]);
                }
            }
            None => out = Some(a[i]),
        }
        i += 1;
    }
    out
}

/// Gets the element-wise maximum between `a` and `b`.
pub(crate) fn elem_max<const N: usize, T: AnyInt>(a: &[T; N], b: &[T; N]) -> [T; N] {
    let mut out = [T::default(); N];
    for i in 0..N {
        out[i] = a[i].max(b[i]);
    }
    out
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_make_uppercase() {
        assert_eq!(make_uppercase(b""), *b"");
        assert_eq!(make_uppercase(b"A"), *b"A");
        assert_eq!(make_uppercase(b"a"), *b"A");
        assert_eq!(make_uppercase(b"1"), *b"1");
        assert_eq!(make_uppercase(b"ACGTN1acgtn1 "), *b"ACGTN1ACGTN1 ");
    }

    #[test]
    fn test_position() {
        assert_eq!(position(b"", b'A'), None);
        assert_eq!(position(b"A", b'A'), Some(0));
        assert_eq!(position(b"a", b'A'), None);
        assert_eq!(position(b"Aa1", b'a'), Some(1));
        assert_eq!(position(b"ACGTN1acgtn1 ", b'n'), Some(10));
    }

    #[test]
    fn test_is_unique() {
        assert!(is_unique(b""), "Should be true!");
        assert!(is_unique(b"A"), "Should be true!");
        assert!(!is_unique(b"AA"), "Should be false!");
        assert!(is_unique(b"ACGTN"), "Should be true!");
        assert!(!is_unique(b"ACGTNA"), "Should be false!");
    }

    #[test]
    fn test_is_subset_with_subset_present() {
        let arr: [u8; 5] = [1, 2, 3, 4, 5];
        let subset: [u8; 3] = [4, 2, 5];
        assert!(is_subset(&arr, &subset));
    }

    #[test]
    fn test_is_subset_with_subset_not_present() {
        let arr: [u8; 5] = [1, 2, 3, 4, 5];
        let subset: [u8; 3] = [6, 7, 8];
        assert!(!is_subset(&arr, &subset));
    }

    #[test]
    fn test_is_subset_with_empty_subset() {
        let arr: [u8; 5] = [1, 2, 3, 4, 5];
        let subset: [u8; 0] = [];
        assert!(is_subset(&arr, &subset));
    }

    #[test]
    fn test_is_subset_with_empty_array() {
        let arr: [u8; 0] = [];
        let subset: [u8; 0] = [];
        assert!(is_subset(&arr, &subset));

        let subset: [u8; 1] = [1];
        assert!(!is_subset(&arr, &subset));
    }

    #[test]
    fn test_is_subset_with_duplicate_elements_in_subset() {
        let arr: [u8; 5] = [1, 2, 2, 3, 4];
        let subset: [u8; 3] = [2, 2, 3];
        assert!(is_subset(&arr, &subset));
    }

    #[test]
    fn test_is_subset_with_all_elements_present() {
        let arr: [u8; 5] = [1, 2, 3, 4, 5];
        let subset: [u8; 5] = [2, 3, 1, 4, 5];
        assert!(is_subset(&arr, &subset));
    }

    #[test]
    fn test_is_subset_with_some_elements_missing() {
        let arr: [u8; 5] = [1, 2, 3, 4, 5];
        let subset: [u8; 4] = [1, 2, 6, 4];
        assert!(!is_subset(&arr, &subset));
    }

    #[test]
    fn test_elem_max_positive_integers() {
        let a = [1, 3, 5, 7];
        let b = [2, 2, 6, 4];
        let result = elem_max(&a, &b);
        assert_eq!(result, [2, 3, 6, 7]);
    }

    #[test]
    fn test_elem_max_negative_integers() {
        let a = [-1, -3, -5, -7];
        let b = [-2, -2, -6, -4];
        let result = elem_max(&a, &b);
        assert_eq!(result, [-1, -2, -5, -4]);
    }

    #[test]
    fn test_elem_max_mixed_integers() {
        let a = [-1, 3, -5, 7];
        let b = [2, -2, 6, -4];
        let result = elem_max(&a, &b);
        assert_eq!(result, [2, 3, 6, 7]);
    }

    #[test]
    fn test_elem_max_equal_arrays() {
        let a = [5, 5, 5, 5];
        let b = [5, 5, 5, 5];
        let result = elem_max(&a, &b);
        assert_eq!(result, [5, 5, 5, 5]);
    }

    #[test]
    fn test_arr_max_typical_values() {
        let a = [1, 3, 5, 2, 4];
        let result = arr_max(&a);
        assert_eq!(result, Some(5));
    }

    #[test]
    fn test_arr_max_identical_values() {
        let a = [7, 7, 7, 7];
        let result = arr_max(&a);
        assert_eq!(result, Some(7));
    }

    #[test]
    fn test_arr_max_single_element() {
        let a = [42];
        let result = arr_max(&a);
        assert_eq!(result, Some(42));
    }

    #[test]
    fn test_arr_max_with_max_u8() {
        let a = [0, 255, 128, 64];
        let result = arr_max(&a);
        assert_eq!(result, Some(255));
    }

    #[test]
    fn test_arr_max_empty_array() {
        let a: &[u8; 0] = &[];
        let result = arr_max(a);
        assert_eq!(result, None);
    }
}
