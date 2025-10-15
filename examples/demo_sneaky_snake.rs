use zoe::alignment::sneaky_snake;

fn main() {
    let reference: &[u8] = b"GGTGCAGAGCTC";
    let query: &[u8] = b"GGTGAGAGTTGT";
    let threshold: f32 = 3. / 12.;

    match sneaky_snake(reference, query, threshold) {
        Some(true) => println!("The sequences are similar within the given threshold."),
        Some(false) => println!("The sequences are not similar within the given threshold."),
        None => println!("Invalid input or threshold."),
    }
}
