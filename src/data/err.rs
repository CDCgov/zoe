use std::fmt::Display;

#[macro_export]
macro_rules! unwrap_or_return_some_err {
    ($expression:expr) => {
        match $expression {
            Ok(v) => v,
            Err(e) => return Some(Err(e)),
        }
    };
}

/// Trait for specifying getting exit codes from errors.
pub trait GetCode {
    fn get_code(&self) -> i32 {
        1
    }
}

impl GetCode for std::io::Error {
    #[must_use]
    #[inline]
    fn get_code(&self) -> i32 {
        self.raw_os_error().unwrap_or(1)
    }
}

/// Trait for providing more graceful [`expect()`](std::result::Result::expect)
/// behavior but with a status code provided by [`GetCode`].
pub trait OrFail<T> {
    fn unwrap_or_fail(self) -> T;
    fn unwrap_or_die(self, msg: &str) -> T;
}

impl<T, E> OrFail<T> for Result<T, E>
where
    E: GetCode + Display,
{
    fn unwrap_or_fail(self) -> T {
        match self {
            Ok(result) => result,
            Err(e) => {
                eprintln!("Error: {e}");
                std::process::exit(e.get_code());
            }
        }
    }

    fn unwrap_or_die(self, msg: &str) -> T {
        match self {
            Ok(result) => result,
            Err(e) => {
                eprintln!("Error: {msg}\n\n{e}");
                std::process::exit(e.get_code());
            }
        }
    }
}

#[derive(Debug)]
pub enum DistanceError {
    NoData,
    NotComparable,
}

impl Display for DistanceError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let message = match self {
            DistanceError::NoData => "No data found in at least one function argument",
            DistanceError::NotComparable => "Data only had invalid state that was not able to be compared",
        };

        write!(f, "{message}")
    }
}

impl std::error::Error for DistanceError {}

#[derive(Debug)]
pub enum AlignmentError {
    EmptyQuery,
    BadGapWeights,
}

impl Display for AlignmentError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let message = match self {
            AlignmentError::EmptyQuery => "The alignment query was empty",
            AlignmentError::BadGapWeights => "The gap extend penalty was greater than the gap open penalty",
        };

        write!(f, "{message}")
    }
}

impl std::error::Error for AlignmentError {}
