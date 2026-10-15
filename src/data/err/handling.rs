//! A module containing traits, macros, and functions for handling errors or
//! results when they occur.

use crate::data::err::{DisplayErrStack, ErrorWithContext, MAX_ERROR_CHAIN_DEPTH};
use std::{error::Error, hint::cold_path};

/// A trait for providing more graceful error reporting and aborting. For
/// similar methods on [`Result`], see [`OrFail`].
///
/// A raw OS error code is used when one is available in the source chain;
/// otherwise, the process exits with code `1`. Any context available in
/// [`Error::source`] is displayed.
pub trait Fail {
    /// Exits the program, writing the error and any information in
    /// [`Error::source`] to stderr.
    fn fail(self) -> !;

    /// Exits the program, writing the provided message, the error, and any
    /// information in [`Error::source`] to stderr.
    fn die(self, msg: &str) -> !;
}

impl<E> Fail for E
where
    E: Error + 'static,
{
    #[cold]
    fn fail(self) -> ! {
        let code = io_exit_code(&self);
        exit_with_error(&self, None, code)
    }

    #[cold]
    fn die(self, msg: &str) -> ! {
        let code = io_exit_code(&self);
        exit_with_error(&self, Some(msg), code)
    }
}

/// A trait for providing more graceful error reporting and aborting.
pub trait OrFail<T> {
    /// Unwraps the result, writing the error and any information in
    /// [`Error::source`] to stderr.
    ///
    /// A raw OS error code is used as the exit code when one is available in
    /// the source chain; otherwise, the process exits with code `1`.
    fn unwrap_or_fail(self) -> T;

    /// Unwraps the result, writing the provided message, the error, and any
    /// information in [`Error::source`] to stderr.
    ///
    /// A raw OS error code is used as the exit code when one is available in
    /// the source chain; otherwise, the process exits with code `1`.
    fn unwrap_or_die(self, msg: &str) -> T;

    /// Unwraps the result, writing the error and any information in
    /// [`Error::source`] to stderr and exiting with `code` on failure.
    fn unwrap_or_exit(self, code: i32) -> T;
}

impl<T, E> OrFail<T> for Result<T, E>
where
    E: Error + 'static,
{
    fn unwrap_or_fail(self) -> T {
        match self {
            Ok(result) => result,
            Err(e) => {
                cold_path();
                e.fail()
            }
        }
    }

    fn unwrap_or_die(self, msg: &str) -> T {
        match self {
            Ok(result) => result,
            Err(e) => {
                cold_path();
                e.die(msg)
            }
        }
    }

    fn unwrap_or_exit(self, code: i32) -> T {
        match self {
            Ok(result) => result,
            Err(error) => {
                cold_path();
                exit_with_error(&error, None, code)
            }
        }
    }
}

/// Finds the first raw OS error code in an error's chain.
///
/// Unlike [`Error::source`], this descends into transparent sources and IO
/// payloads. These can hold a raw OS error code which is not visible by the
/// wrapping error.
fn io_exit_code(error: &(dyn Error + 'static)) -> i32 {
    let mut current = Some(error);

    for _ in 0..MAX_ERROR_CHAIN_DEPTH {
        let Some(error) = current else {
            break;
        };

        // Extract the raw_os_error if present
        if let Some(io_error) = error.downcast_ref::<std::io::Error>()
            && let Some(code) = io_error.raw_os_error()
        {
            return code;
        }

        current = if let Some(error) = error.downcast_ref::<ErrorWithContext>()
            && let Some(error) = error.transparent_err()
        {
            // Get transparently wrapped error
            Some(error)
        } else if let Some(payload) = error.downcast_ref::<std::io::Error>().and_then(std::io::Error::get_ref) {
            // Get IO payload
            Some(payload)
        } else {
            // Get source error
            error.source()
        };
    }

    1
}

/// Writes an error stack to stderr and exits with `code`.
fn exit_with_error(error: &(dyn Error + 'static), msg: Option<&str>, code: i32) -> ! {
    match (std::env::current_exe(), msg) {
        (Ok(bin), Some(msg)) => eprintln!("Error in {bin}: {msg}", bin = bin.display()),
        (Ok(bin), None) => eprintln!("Error in {bin}", bin = bin.display()),
        (Err(_), Some(msg)) => eprintln!("Error: {msg}"),
        (Err(_), None) => eprintln!("Error in program"),
    }

    eprint!("{}", error.display_stack());
    std::process::exit(code);
}

/// A macro for unwrapping a [`Result`] and propagating any error as a
/// `Some(Err(e))`.
///
/// This is especially useful for fallible iterators, where results need to be
/// wrapped in [`Some`].
#[macro_export]
macro_rules! unwrap_or_return_some_err {
    ($expression:expr) => {
        match $expression {
            Ok(v) => v,
            Err(e) => return Some(Err(e)),
        }
    };
}
