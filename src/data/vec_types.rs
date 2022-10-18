#[inline]
pub(crate) fn find_and_replace<T>(v: &mut [T], needle: T, replacement: T)
where
    T: Copy + PartialEq,
{
    for b in v.iter_mut() {
        if *b == needle {
            *b = replacement;
        }
    }
}
