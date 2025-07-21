use std::error::Error;
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
    /// Retrieves the exit code associated with a given error.
    ///
    /// The blanket implementation returns `1`. We also implement on
    /// [`std::io::Error`] to return the [`raw_os_error`] if available.
    ///
    /// [`raw_os_error`]: std::io::Error::raw_os_error
    #[inline]
    #[must_use]
    fn get_code(&self) -> i32 {
        1
    }
}

impl GetCode for std::io::Error {
    #[inline]
    fn get_code(&self) -> i32 {
        if let Some(code) = self.raw_os_error() {
            return code;
        }

        let mut source = self.source();
        while let Some(err) = source {
            if let Some(e) = err.downcast_ref::<std::io::Error>()
                && let Some(code) = e.raw_os_error()
            {
                return code;
            }
            source = err.source();
        }

        1
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
pub enum QueryProfileError {
    EmptyQuery,
    GapOpenOutOfRange { gap_open: i8 },
    GapExtendOutOfRange { gap_extend: i8 },
    BadGapWeights { gap_open: i8, gap_extend: i8 },
}

impl Display for QueryProfileError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let message = match self {
            QueryProfileError::EmptyQuery => "The alignment query was empty",
            QueryProfileError::GapOpenOutOfRange { gap_open } => {
                &format!("The gap open weight must be between -127 and 0, but {gap_open} was provided")
            }
            QueryProfileError::GapExtendOutOfRange { gap_extend } => {
                &format!("The gap extend weight must be between -127 and 0, but {gap_extend} was provided")
            }
            QueryProfileError::BadGapWeights { gap_open, gap_extend } => &format!(
                "The gap open weight must be less than or equal to the gap extend weight, but {gap_open} (gap open) and {gap_extend} (gap extend) were provided"
            ),
        };

        write!(f, "{message}")
    }
}

impl std::error::Error for QueryProfileError {}
