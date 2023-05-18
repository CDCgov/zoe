/// Takes square root using the Babylonian Method using 40 iterations
// Relevant discussion: <https://www.codeproject.com/Articles/69941/Best-Square-Root-Method-Algorithm-Function-Precisi>
pub(crate) const fn sqrt_baby(n: f64) -> f64 {
    let mut i: f64 = 0.0;
    while i * i <= n {
        i += 0.1_f64;
    }

    let mut x1 = i;
    let mut x2 = 0.0;

    let mut j: usize = 0;
    while j < 40 {
        x2 = ((n / x1) + x1) / 2.0;
        x1 = x2;
        j += 1;
    }
    x2
}
