// TODO: Make into a trait when functionality becomes available.
// See: https://github.com/rust-lang/rust-project-goals/issues/106

pub(crate) const fn make_uppercase<const N: usize>(a: &[u8; N]) -> [u8; N] {
    let mut b = [0; N];
    let mut i = 0;
    while i < a.len() {
        b[i] = a[i].to_ascii_uppercase();
        i += 1;
    }
    b
}

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
    fn test_is_uniqe() {
        assert!(is_unique(b""), "Should be true!");
        assert!(is_unique(b"A"), "Should be true!");
        assert!(!is_unique(b"AA"), "Should be false!");
        assert!(is_unique(b"ACGTN"), "Should be true!");
        assert!(!is_unique(b"ACGTNA"), "Should be false!");
    }
}
