// Make into a trait when functionality becomes available.
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_is_uniqe() {
        assert!(is_unique(b"ACGTN"), "Should be true!");
        assert!(!is_unique(b"ACGTNA"), "Should be false!");
    }
}
