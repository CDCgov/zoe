//! Error handling utilities for *Zoe*.
//!
//! This module provides:
//!
//! - [`ErrorWithContext`] and [`WithErrorContext`] for wrapping errors with
//!   additional context while preserving the error source chain.
//! - [`GetCode`] and [`OrFail`] for graceful CLI error handling with exit
//!   codes.
//! - [`open_nonempty_file`] for opening files with automatic path context on
//!   errors.
//!
//! ## Error Handling Philosophy
//!
//! *Zoe* often wraps [`std::io::Error`] for I/O and parsing operations rather
//! than defining new error types. Domain-specific errors (e.g., alignment,
//! distance) use enums in their respective modules. [`std::fmt::Display`] is
//! implemented only at the immediate error level—callers and loggers must
//! iterate through the source chain via [`std::error::Error::source`], use
//! Zoe's [`OrFail`], or an external crates like `anyhow`.
//!
//! [`ErrorWithContext`]: crate::data::err::ErrorWithContext
//! [`WithErrorContext`]: crate::data::err::WithErrorContext
//! [`GetCode`]: crate::data::err::GetCode
//! [`OrFail`]: crate::data::err::OrFail
//! [`open_nonempty_file`]: crate::data::err::open_nonempty_file

use std::{
    error::Error,
    fmt::{Display, Write},
    fs::File,
    io::ErrorKind,
    path::Path,
};

#[macro_export]
macro_rules! unwrap_or_return_some_err {
    ($expression:expr) => {
        match $expression {
            Ok(v) => v,
            Err(e) => return Some(Err(e)),
        }
    };
}

/// Trait for specifying getting exit codes originating from IO errors.
///
/// Implementing this trait allows the error type to work with [`OrFail`]. If
/// this is being implemented for a top-level error containing a nested error,
/// one should manually implement [`get_code`] to retrieve the underlying code
/// for the nested error, rather than using the blanket implementation.
///
/// [`get_code`]: GetCode::get_code
pub trait GetCode {
    /// Retrieves the exit code associated with a given error.
    ///
    /// ## Validity
    ///
    /// If this method is manually implemented, then it must recursively call
    /// [`get_code`] on [`Error::source`]. Any other behavior or logic is not
    /// guaranteed to be consistent within *Zoe*, and may or may not be
    /// correctly applied when using wrapped errors ([`ErrorWithContext`]).
    ///
    /// The blanket implementation returns `1`. We also implement on
    /// [`std::io::Error`] to return the [`raw_os_error`] if available.
    ///
    /// [`raw_os_error`]: std::io::Error::raw_os_error
    /// [`get_code`]: GetCode::get_code
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

/// A trait for providing more graceful error reporting and aborting.
///
/// A status code is provided by [`GetCode`], and any context available in
/// [`Error::source`] is displayed.
///
/// <div class="warning note">
///
/// **Note**
///
/// To get full utility out of this trait, custom top-level errors should
/// manually implement [`std::error::Error::source`] and get whatever field or
/// variants contains the nested errors. In addition, [`GetCode`] should
/// likewise be implemented manually to retrieve the underlying codes for nested
/// [`std::io::Error`].
///
///
/// </div>
pub trait OrFail<T> {
    /// Unwraps the result, writing the error and any information in
    /// [`Error::source`] to stderr.
    fn unwrap_or_fail(self) -> T;

    /// Unwraps the result, writing the provided message, the error, and any
    /// information in [`Error::source`] to stderr.
    fn unwrap_or_die(self, msg: &str) -> T;
}

impl<T, E> OrFail<T> for Result<T, E>
where
    E: GetCode + Display + Error,
{
    fn unwrap_or_fail(self) -> T {
        match self {
            Ok(result) => result,
            Err(e) => {
                if let Ok(bin) = std::env::current_exe() {
                    eprintln!("Error in {b}", b = bin.display());
                    eprintln!("  → {e}", e = IndentWrapper(&e));
                } else {
                    eprintln!("Error: {e}");
                }
                let mut source = e.source();
                while let Some(err) = source {
                    eprintln!("  → {err}", err = IndentWrapper(err));
                    source = err.source();
                }
                std::process::exit(e.get_code());
            }
        }
    }

    fn unwrap_or_die(self, msg: &str) -> T {
        match self {
            Ok(result) => result,
            Err(e) => {
                if let Ok(bin) = std::env::current_exe() {
                    eprintln!("Error in {b}: {msg}\n\n{e}", b = bin.display());
                } else {
                    eprintln!("Error: {msg}\n\n{e}");
                }
                let mut source = e.source();
                while let Some(err) = source {
                    eprintln!("  → {err}", err = IndentWrapper(err));
                    source = err.source();
                }
                std::process::exit(e.get_code());
            }
        }
    }
}

/// A wrapper around an error with a new message, and the original error
/// accessible via [`Error::source`].
///
/// When [`unwrap_or_fail`] or [`unwrap_or_die`] is used, this will display the
/// information for the original error as well.
///
/// [`unwrap_or_fail`]: OrFail::unwrap_or_fail
/// [`unwrap_or_die`]: OrFail::unwrap_or_die
#[must_use]
#[derive(Debug)]
pub struct ErrorWithContext {
    description: String,
    source:      Box<dyn Error + Send + Sync>,
}

impl std::fmt::Display for ErrorWithContext {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.description)
    }
}

impl Error for ErrorWithContext {
    #[inline]
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(self.source.as_ref())
    }
}

impl From<ErrorWithContext> for std::io::Error {
    #[inline]
    fn from(e: ErrorWithContext) -> Self {
        std::io::Error::other(e)
    }
}

/// An extension trait for [`Error`] allowing additional context to be added via
/// a [`ErrorWithContext`].
pub trait WithErrorContext {
    /// Wraps the error in an [`ErrorWithContext`] with the given description.
    ///
    /// The `description` field may be anything implementing `Into<String>`.
    /// Passing an owned `String` avoids an extra allocation.
    fn with_context(self, description: impl Into<String>) -> ErrorWithContext;

    /// Convenience function for adding type context.
    fn with_type_context<T>(self) -> ErrorWithContext;

    /// Convenience function for adding file context.
    ///
    /// The context will be formatted as `msg: file`. The `msg` field may be
    /// anything implementing [`Display`].
    fn with_file_context(self, msg: impl Display, file: impl AsRef<Path>) -> ErrorWithContext;
}

impl<E: Error + Send + Sync + 'static> WithErrorContext for E {
    #[inline]
    fn with_context(self, description: impl Into<String>) -> ErrorWithContext {
        ErrorWithContext {
            description: description.into(),
            source:      Box::new(self),
        }
    }

    #[inline]
    fn with_type_context<T>(self) -> ErrorWithContext {
        let name = std::any::type_name::<T>();
        let description = format!(
            "Failure in {}",
            name.split('<').next().unwrap_or(name).rsplit("::").next().unwrap_or(name)
        );

        ErrorWithContext {
            description,
            source: Box::new(self),
        }
    }

    #[inline]
    fn with_file_context(self, msg: impl Display, file: impl AsRef<Path>) -> ErrorWithContext {
        Self::with_context(self, format!("{msg}: '{path}'", path = file.as_ref().display()))
    }
}

/// An extension trait for [`Result`] allowing additional context to be added to
/// an [`Err`] variant via a [`ErrorWithContext`].
///
/// The methods are similar to [`WithErrorContext`], but are implemented for
/// results.
pub trait ResultWithErrorContext {
    /// The type of the [`Ok`] variant in the result.
    type Ok;

    /// Wraps the [`Err`] variant in a [`ErrorWithContext`] with the given
    /// description.
    ///
    /// ## Errors
    ///
    /// Propagates errors in `self`, with the added context.
    fn with_context(self, description: impl Into<String>) -> Result<Self::Ok, ErrorWithContext>;

    /// Wraps the [`Err`] variant in a [`ErrorWithContext`] by adding type
    /// context.
    ///
    /// ## Errors
    ///
    /// Propagates errors in `self`, with the added context.
    fn with_type_context<T>(self) -> Result<Self::Ok, ErrorWithContext>;

    /// Wraps the [`Err`] variant in a [`ErrorWithContext`] by adding file
    /// context.
    ///
    /// ## Errors
    ///
    /// Propagates errors in `self`, with the added context.
    fn with_file_context(self, msg: impl Display, file: impl AsRef<Path>) -> Result<Self::Ok, ErrorWithContext>;
}

impl<Ok, E: Error + Send + Sync + 'static> ResultWithErrorContext for Result<Ok, E> {
    type Ok = Ok;

    #[inline]
    fn with_context(self, description: impl Into<String>) -> Result<Ok, ErrorWithContext> {
        self.map_err(|e| e.with_context(description))
    }

    #[inline]
    fn with_type_context<T>(self) -> Result<Ok, ErrorWithContext> {
        self.map_err(WithErrorContext::with_type_context::<T>)
    }

    #[inline]
    fn with_file_context(self, msg: impl Display, file: impl AsRef<Path>) -> Result<Ok, ErrorWithContext> {
        self.map_err(|e| e.with_file_context(msg, file))
    }
}

/// Opens a file, checking that it is non-empty, and wrapping any errors with
/// the file path for context.
///
/// ## Errors
///
/// Returns an error if:
///
/// - The file cannot be opened
/// - File metadata cannot be read
/// - The file is empty
///
/// All errors include the file path in the error message and preserve the
/// original error via [`Error::source`].
#[inline]
pub fn open_nonempty_file(path: impl AsRef<Path>) -> std::io::Result<File> {
    let path = path.as_ref();

    let file = File::open(path)
        .map_err(|e| std::io::Error::other(e.with_context(format!("Failed to open file: '{}'", path.display()))))?;

    let metadata = file
        .metadata()
        .map_err(|e| std::io::Error::other(e.with_context(format!("Failed to read metadata: '{}'", path.display()))))?;

    if metadata.len() == 0 {
        return Err(std::io::Error::new(
            ErrorKind::InvalidInput,
            format!("File is empty: '{}'", path.display()),
        ));
    }

    Ok(file)
}

/// A wrapper around [`std::fmt::Formatter`] which automatically indents all new
/// lines by four spaces.
///
/// This is a helper struct for the formatting used in [`unwrap_or_fail`] and
/// [`unwrap_or_die`].
///
/// [`unwrap_or_fail`]: OrFail::unwrap_or_fail
/// [`unwrap_or_die`]: OrFail::unwrap_or_die
struct IndentFormatter<'a, 'b>(&'a mut std::fmt::Formatter<'b>);

impl Write for IndentFormatter<'_, '_> {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
        let mut parts = s.split('\n');
        let Some(first_part) = parts.next() else { return Ok(()) };
        self.0.write_str(first_part)?;

        for part in parts {
            self.0.write_str("\n    ")?;
            self.0.write_str(part)?;
        }

        Ok(())
    }

    fn write_char(&mut self, c: char) -> std::fmt::Result {
        if c == '\n' {
            self.0.write_str("\n    ")
        } else {
            self.0.write_char(c)
        }
    }
}

/// A wrapper type altering the implementation of [`Display`], such that any new
/// lines are automatically indented with four spaces.
///
/// This is a helper struct for the formatting used in [`unwrap_or_fail`] and
/// [`unwrap_or_die`].
///
/// [`unwrap_or_fail`]: OrFail::unwrap_or_fail
/// [`unwrap_or_die`]: OrFail::unwrap_or_die
struct IndentWrapper<T>(T);

impl<T: Display> Display for IndentWrapper<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(IndentFormatter(f), "{}", self.0)
    }
}
