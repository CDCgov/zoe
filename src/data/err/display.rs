//! A module with code for displaying errors and error stacks.

use crate::data::err::MAX_ERROR_CHAIN_DEPTH;
use std::{
    error::Error,
    fmt::{Display, Write},
};

/// An extension trait for an [`Error`] enabling it to be displayed alongside
/// its error stack using [`Error::source`].
///
/// If aborting, consider using [`fail`] or [`die`]. Otherwise, this can be
/// helpful in displaying warnings.
///
/// [`fail`]: crate::data::err::Fail::fail
/// [`die`]: crate::data::err::Fail::die
pub trait DisplayErrStack {
    /// Returns a displayable representation of the error alongside its error
    /// stack using [`Error::source`].
    ///
    /// This using `→` and two spaces of indent before each item, and includes a
    /// newline at the end.
    fn display_stack(&self) -> ErrStackDisplay<'_>;
}

impl<E> DisplayErrStack for E
where
    E: Error + 'static,
{
    fn display_stack(&self) -> ErrStackDisplay<'_> {
        ErrStackDisplay(self)
    }
}

impl DisplayErrStack for dyn Error + 'static {
    fn display_stack(&self) -> ErrStackDisplay<'_> {
        ErrStackDisplay(self)
    }
}

impl DisplayErrStack for dyn Error + Send + 'static {
    fn display_stack(&self) -> ErrStackDisplay<'_> {
        ErrStackDisplay(self)
    }
}

impl DisplayErrStack for dyn Error + Sync + 'static {
    fn display_stack(&self) -> ErrStackDisplay<'_> {
        ErrStackDisplay(self)
    }
}

impl DisplayErrStack for dyn Error + Send + Sync + 'static {
    fn display_stack(&self) -> ErrStackDisplay<'_> {
        ErrStackDisplay(self)
    }
}

/// A display wrapper around an error that shows the error and its sources (with
/// [`Error::source`]) in a list, using `→` and two spaces of indent before each
/// item.
///
/// This includes a newline at the end.
pub struct ErrStackDisplay<'a>(&'a (dyn Error + 'static));

impl Display for ErrStackDisplay<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut maybe_err = Some(self.0);

        for _ in 0..MAX_ERROR_CHAIN_DEPTH {
            let Some(err) = maybe_err else {
                return Ok(());
            };

            writeln!(
                f,
                "  → {err}",
                err = IndentWrapper {
                    val:    err,
                    indent: "    ",
                }
            )?;

            maybe_err = err.source();
        }

        if maybe_err.is_some() {
            writeln!(f, "  → [error source chain truncated]")?;
        }

        Ok(())
    }
}

/// A wrapper type altering the implementation of [`Display`], such that any new
/// lines are automatically indented with a specified string.
///
/// This is a helper struct for the formatting used in [`fail`] and [`die`].
///
/// [`fail`]: crate::data::err::Fail::fail
/// [`die`]: crate::data::err::Fail::die
pub(super) struct IndentWrapper<T> {
    /// The value to display.
    pub val:    T,
    /// The string to use when indenting lines after the first.
    pub indent: &'static str,
}

impl<T: Display> Display for IndentWrapper<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            IndentFormatter {
                formatter: f,
                indent:    self.indent,
            },
            "{}",
            self.val
        )
    }
}

/// A wrapper around [`std::fmt::Formatter`] which automatically indents all new
/// lines with a specified string.
///
/// This is a helper struct for the formatting used in [`fail`] and [`die`].
///
/// [`fail`]: crate::data::err::Fail::fail
/// [`die`]: crate::data::err::Fail::die
struct IndentFormatter<'a, 'b> {
    formatter: &'a mut std::fmt::Formatter<'b>,
    indent:    &'static str,
}

impl Write for IndentFormatter<'_, '_> {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
        let mut parts = s.split('\n');
        let Some(first_part) = parts.next() else { return Ok(()) };
        self.formatter.write_str(first_part)?;

        for part in parts {
            self.formatter.write_char('\n')?;
            self.formatter.write_str(self.indent)?;
            self.formatter.write_str(part)?;
        }

        Ok(())
    }

    fn write_char(&mut self, c: char) -> std::fmt::Result {
        if c == '\n' {
            self.formatter.write_char('\n')?;
            self.formatter.write_str(self.indent)
        } else {
            self.formatter.write_char(c)
        }
    }
}
