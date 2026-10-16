//! A module defining [`ErrorWithContext`] and related constructors/traits.

use crate::data::err::{DisplayErrStack, IndentWrapper};
use std::{
    error::Error,
    fmt::{Debug, Display},
    hint::cold_path,
    path::Path,
};

/// An error type supporting context and a backtrace.
///
/// Specifically, this error can hold up to three things:
///
/// 1. A line of context describing the error. This appears as one item in the
///    [`OrFail`] backtrace.
/// 2. Any subitems (additional indented lines with more information that appear
///    below the line of context). This is useful for including the values of
///    variables or other useful information.
/// 3. An optional source error, which this error wraps. Using
///    [`unwrap_or_fail`] or [`unwrap_or_die`] cause the source error to be
///    shown in the backtrace. This source is accessible via [`Error::source`].
///
/// This can be converted to [`std::io::Error`] with [`Into`]. Hence, in
/// functions returning [`std::io::Result`], the `?` operator can be used after
/// adding context.
///
/// A [`std::io::Error`] can also be converted into a [`ErrorWithContext`]
/// without adding additional context via [`From`]. The implementation ensures
/// that the error does not appear twice in the public source chain. This can be
/// useful to return `Result<(), ErrorWithContext>` from the `main` function,
/// which will automatically format any errors and the source chain upon
/// failure.
///
/// [`OrFail`]: crate::data::err::OrFail
/// [`with_subitem`]: WithSubitem::with_subitem
/// [`unwrap_or_fail`]: crate::data::err::OrFail::unwrap_or_fail
/// [`unwrap_or_die`]: crate::data::err::OrFail::unwrap_or_die
#[must_use]
pub struct ErrorWithContext {
    /// The inner representation. Using a single fat pointer is better than
    /// storing the description and source fields directly, since it minimizes
    /// the size of [`ErrorWithContext`] and hence the size of `Result<T,
    /// ErrorWithContext>`.
    repr: Box<ErrorRepr>,
}

impl ErrorWithContext {
    /// Constructs a new [`ErrorWithContext`] with the given description,
    /// without a source error or any subitems.
    ///
    /// The `description` may be anything implementing `Into<String>`. Passing
    /// an owned `String` avoids an extra allocation.
    pub fn new(description: impl Into<String>) -> Self {
        ErrorWithContext {
            repr: Box::new(ErrorRepr::WithContext {
                description: description.into(),
                subitem:     None,
                source:      None,
            }),
        }
    }

    /// Constructs a new [`ErrorWithContext`] with the given description and
    /// path, without a source error or any subitems.
    ///
    /// The `description` may be anything implementing `Display`. An allocation
    /// will occur for the resulting error message.
    pub fn new_with_path(description: impl Display, path: impl AsRef<Path>) -> Self {
        let description = format!("{description}: '{path}'", path = path.as_ref().display());

        ErrorWithContext {
            repr: Box::new(ErrorRepr::WithContext {
                description,
                subitem: None,
                source: None,
            }),
        }
    }

    /// If the [`ErrorWithContext`] holds wraps an error without adding an
    /// additional line of context, then the error is returned.
    #[inline]
    #[must_use]
    pub(crate) fn transparent_err(&self) -> Option<&(dyn Error + Send + Sync + 'static)> {
        match self.repr.as_ref() {
            ErrorRepr::WithContext { .. } => None,
            ErrorRepr::NoContext { error, .. } => Some(&**error),
        }
    }
}

impl Debug for ErrorWithContext {
    /// Formats the value using the debug formatter.
    ///
    /// By default, this will print the binary name (if available) followed by a
    /// formatted backtrace of the error. If using the alternate display with
    /// `{:#}`, a traditional debug format is used.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if f.alternate() {
            f.debug_struct("ErrorWithContext").field("repr", &self.repr).finish()
        } else {
            if let Ok(bin) = std::env::current_exe() {
                writeln!(f, "A problem occurred in {b}", b = bin.display())?;
            } else {
                writeln!(f)?;
            }

            write!(f, "{}", self.display_stack())
        }
    }
}

/// The inner representation for an [`ErrorWithContext`]. This is either a line
/// of context with optional subitems and an optional source, or a wrapper
/// without additional context around a different error.
///
/// This is wrapped in a [`Box`] in [`ErrorWithContext`] to reduce the memory of
/// `Result<T, ErrorWithContext>` in the `Ok` case.
#[derive(Debug)]
enum ErrorRepr {
    /// An error including a line of context, along with possibly any subitems
    /// or a source error.
    WithContext {
        description: String,
        subitem:     Option<String>,
        source:      Option<Box<dyn Error + Send + Sync>>,
    },
    /// A transparent wrapper around another error without another line of
    /// context. Subitems may still be added.
    NoContext {
        subitem: Option<String>,
        error:   Box<dyn Error + Send + Sync>,
    },
}

impl Display for ErrorRepr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let subitem = match self {
            ErrorRepr::WithContext {
                description, subitem, ..
            } => {
                write!(f, "{description}")?;
                subitem
            }
            ErrorRepr::NoContext { subitem, error } => {
                write!(f, "{error}")?;
                subitem
            }
        };

        if let Some(subitem) = subitem {
            write!(
                f,
                "\n| {}",
                IndentWrapper {
                    val:    subitem,
                    indent: "| ",
                }
            )?;
        }

        Ok(())
    }
}

impl Display for ErrorWithContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.repr)
    }
}

impl Error for ErrorWithContext {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self.repr.as_ref() {
            ErrorRepr::WithContext {
                source: Some(source), ..
            } => Some(&**source),
            ErrorRepr::WithContext { source: None, .. } => None,
            ErrorRepr::NoContext { error, .. } => error.source(),
        }
    }
}

/// An extension trait for [`Error`] allowing additional context to be added via
/// a [`ErrorWithContext`].
pub trait WithErrorContext {
    /// Wraps the error in an [`ErrorWithContext`] with the given description.
    ///
    /// The `description` may be anything implementing `Into<String>`. Passing
    /// an owned `String` avoids an extra allocation.
    fn with_context(self, description: impl Into<String>) -> ErrorWithContext;

    /// Wraps the error in an [`ErrorWithContext`] by adding type context.
    fn with_type_context<T>(self) -> ErrorWithContext;

    /// Wraps the error in an [`ErrorWithContext`] by adding path context.
    ///
    /// The context will be formatted as `msg: 'path'`. The `msg` may be
    /// anything implementing [`Display`].
    fn with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> ErrorWithContext;
}

impl<E: Error + Send + Sync + 'static> WithErrorContext for E {
    // Do not inline, since this is cold code
    fn with_context(self, description: impl Into<String>) -> ErrorWithContext {
        ErrorWithContext {
            repr: Box::new(ErrorRepr::WithContext {
                description: description.into(),
                subitem:     None,
                source:      Some(Box::new(self)),
            }),
        }
    }

    // Do not inline, since this is cold code
    fn with_type_context<T>(self) -> ErrorWithContext {
        let name = std::any::type_name::<T>();
        let description = format!(
            "Failure in {}",
            name.split('<').next().unwrap_or(name).rsplit("::").next().unwrap_or(name)
        );

        Self::with_context(self, description)
    }

    // Do not inline, since this is cold code
    fn with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> ErrorWithContext {
        Self::with_context(self, format!("{msg}: '{path}'", path = file.as_ref().display()))
    }
}

/// An extension trait for [`ErrorWithContext`] allowing an indented subitem to
/// be added to the error (without adding a new error to the backtrace).
pub trait WithSubitem {
    /// Adds a subitem with the given `message` to the error without adding a
    /// new error to the backtrace.
    ///
    /// [`WithErrorContext::with_context`] creates a new entry in the backtrace
    /// (displayed using `→`), whereas this method adds a message that is
    /// indented beneath the error and indicated using `|`. For example:
    ///
    /// ```text
    /// Error in /path/to/binary
    ///   → Failed to load reads from file: input.fastq
    ///   → Failed to deinterleave records due to mismatching headers
    ///     | Header 1: SIM:1:FCX:1:15:6329:1045 1:N:0:2
    ///     | Header 2: SIM:1:FCX:1:15:2345:1001 2:N:0:2
    ///   → x_pos fields did not agree!
    /// ```
    fn with_subitem(self, message: impl Into<String>) -> ErrorWithContext;
}

impl WithSubitem for ErrorWithContext {
    fn with_subitem(mut self, message: impl Into<String>) -> ErrorWithContext {
        let subitem = match self.repr.as_mut() {
            ErrorRepr::WithContext { subitem, .. } | ErrorRepr::NoContext { subitem, .. } => subitem,
        };
        let message = message.into();
        if let Some(subitem) = subitem {
            subitem.push('\n');
            subitem.push_str(&message);
        } else {
            *subitem = Some(message);
        }
        self
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

    /// Wraps the [`Err`] variant in an [`ErrorWithContext`] with the given
    /// description.
    ///
    /// The `description` may be anything implementing `Into<String>`. Passing
    /// an owned `String` avoids an extra allocation.
    ///
    /// ## Errors
    ///
    /// Propagates errors in `self`, with the added context.
    fn with_context(self, description: impl Into<String>) -> Result<Self::Ok, ErrorWithContext>;

    /// Wraps the [`Err`] variant in an [`ErrorWithContext`] by adding type
    /// context.
    ///
    /// ## Errors
    ///
    /// Propagates errors in `self`, with the added context.
    fn with_type_context<T>(self) -> Result<Self::Ok, ErrorWithContext>;

    /// Wraps the [`Err`] variant in an [`ErrorWithContext`] by adding path
    /// context.
    ///
    /// The context will be formatted as `msg: 'path'`. The `msg` may be
    /// anything implementing [`Display`].
    ///
    /// ## Errors
    ///
    /// Propagates errors in `self`, with the added context.
    fn with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> Result<Self::Ok, ErrorWithContext>;
}

impl<Ok, E: WithErrorContext> ResultWithErrorContext for Result<Ok, E> {
    type Ok = Ok;

    #[inline]
    fn with_context(self, description: impl Into<String>) -> Result<Ok, ErrorWithContext> {
        self.map_err(|e| {
            cold_path();
            e.with_context(description)
        })
    }

    #[inline]
    fn with_type_context<T>(self) -> Result<Ok, ErrorWithContext> {
        self.map_err(|e| {
            cold_path();
            e.with_type_context::<T>()
        })
    }

    #[inline]
    fn with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> Result<Ok, ErrorWithContext> {
        self.map_err(|e| {
            cold_path();
            e.with_path_context(msg, file)
        })
    }
}

/// An extension trait for `Result<T, WithErrorContext>` allowing information to
/// be attached to an [`Err`] variant.
///
/// The methods are similar to [`WithSubitem`], but are implemented for results.
pub trait ResultWithSubitem {
    /// Adds a subitem with the given `message` to an [`Err`] variant without
    /// adding a new error to the backtrace.
    ///
    /// [`ResultWithErrorContext::with_context`] creates a new entry in the
    /// backtrace (displayed using `→`), whereas this method adds a message that
    /// is indented beneath the error and indicated using `|`. For example:
    ///
    /// ```text
    /// Error in /path/to/binary
    ///   → Failed to load reads from file: input.fastq
    ///   → Failed to deinterleave records due to mismatching headers
    ///     | Header 1: SIM:1:FCX:1:15:6329:1045 1:N:0:2
    ///     | Header 2: SIM:1:FCX:1:15:2345:1001 2:N:0:2
    ///   → x_pos fields did not agree!
    /// ```
    #[must_use]
    fn with_subitem(self, message: impl Into<String>) -> Self;
}

impl<T> ResultWithSubitem for Result<T, ErrorWithContext> {
    #[inline]
    fn with_subitem(self, message: impl Into<String>) -> Self {
        self.map_err(|e| {
            cold_path();
            e.with_subitem(message)
        })
    }
}

/// An extension trait for [`Error`] allowing it to be converted to
/// [`ErrorWithContext`] without adding any additional context.
pub trait WrapError: Sized + Error + Sync + Send + 'static {
    /// Wraps the error inside an [`ErrorWithContext`] without adding any
    /// context.
    ///
    /// For [`std::io::Error`], prefer calling [`Into::into`].
    fn wrap(self) -> ErrorWithContext {
        ErrorWithContext {
            repr: Box::new(ErrorRepr::NoContext {
                subitem: None,
                error:   Box::new(self),
            }),
        }
    }
}

impl<E: Error + Sync + Send + 'static> WrapError for E {}

/// An extension trait for `Result<T, E>` allowing it to be converted to
/// `Result<T, ErrorWithContext>` without adding any additional context.
pub trait WrapResult {
    /// The type of the [`Ok`] variant in the result.
    type Ok;

    /// Wraps any error inside an [`ErrorWithContext`] without adding any
    /// context.
    ///
    /// If additional context should be added, use [`with_context`]. If the
    /// error type is [`std::io::Error`], prefer calling `map_err(Into::into)`.
    ///
    /// ## Errors
    ///
    /// An error is returned solely if the input result contained an error.
    ///
    /// [`with_context`]: ResultWithErrorContext::with_context
    fn wrap(self) -> Result<Self::Ok, ErrorWithContext>;
}

impl<Ok, E: WrapError> WrapResult for Result<Ok, E> {
    type Ok = Ok;

    fn wrap(self) -> Result<Ok, ErrorWithContext> {
        self.map_err(|e| {
            cold_path();
            e.wrap()
        })
    }
}

impl From<ErrorWithContext> for std::io::Error {
    #[inline]
    fn from(error: ErrorWithContext) -> Self {
        // If the error has no description or subitems and is wrapping a
        // std::io::Error, extract it rather than adding another wrapper. This
        // allows methods like raw_os_error to return Some
        if let ErrorRepr::NoContext { subitem: None, error } = *error.repr {
            match error.downcast::<std::io::Error>() {
                Ok(error) => *error,
                // Downcast failed, so reconstruct the ErrorWithContext and wrap
                // it
                Err(error) => std::io::Error::other(ErrorWithContext {
                    repr: Box::new(ErrorRepr::NoContext { subitem: None, error }),
                }),
            }
        } else {
            std::io::Error::other(error)
        }
    }
}

impl From<std::io::Error> for ErrorWithContext {
    fn from(error: std::io::Error) -> Self {
        // This downcast avoids two wrappers when doing a round-trip conversion
        // (and hence additional Box allocations)
        let error = match error.downcast::<Self>() {
            Ok(error) => return error,
            Err(error) => error,
        };

        ErrorWithContext {
            repr: Box::new(ErrorRepr::NoContext {
                subitem: None,
                error:   Box::new(error),
            }),
        }
    }
}

impl From<std::convert::Infallible> for ErrorWithContext {
    fn from(error: std::convert::Infallible) -> Self {
        match error {}
    }
}
