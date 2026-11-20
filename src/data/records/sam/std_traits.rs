//! Implementations of standard library traits on [`SamData`] and SAM-related
//! types.

use crate::data::{
    err::ResultWithErrorContext,
    sam::{DisplaySamFields, OptArray, SamData, SamDataView, SamOptField, SamOptValue},
};
use std::{fmt::Display, str::FromStr};

impl Display for SamData {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.display())
    }
}

impl Display for SamDataView<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.display())
    }
}

impl Display for SamOptField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}:{}", self.tag[0] as char, self.tag[1] as char, self.value)
    }
}

impl Display for SamOptValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SamOptValue::Char(val) => write!(f, "A:{val}", val = *val as char),
            SamOptValue::Int(val) => write!(f, "i:{val}"),
            SamOptValue::Float(val) => write!(f, "f:{val}"),
            SamOptValue::String(val) => write!(f, "Z:{val}"),
            SamOptValue::Hex(val) => write!(f, "H:{val}"),
            SamOptValue::Array(val) => match val {
                OptArray::I8(vec) => OptArray::fmt_opt_array(f, 'c', vec),
                OptArray::U8(vec) => OptArray::fmt_opt_array(f, 'C', vec),
                OptArray::I16(vec) => OptArray::fmt_opt_array(f, 's', vec),
                OptArray::U16(vec) => OptArray::fmt_opt_array(f, 'S', vec),
                OptArray::I32(vec) => OptArray::fmt_opt_array(f, 'i', vec),
                OptArray::U32(vec) => OptArray::fmt_opt_array(f, 'I', vec),
                OptArray::F32(vec) => OptArray::fmt_opt_array(f, 'f', vec),
            },
        }
    }
}

impl OptArray {
    /// A helper function for displaying a [`SamOptValue::Array`].
    fn fmt_opt_array<T: Display>(f: &mut std::fmt::Formatter<'_>, arr_type: char, vals: &[T]) -> std::fmt::Result {
        write!(f, "B:{arr_type}")?;

        for v in vals {
            write!(f, ",{v}")?;
        }

        Ok(())
    }
}

impl FromStr for SamOptField {
    type Err = std::io::Error;

    fn from_str(field: &str) -> Result<Self, Self::Err> {
        let inv_opt_err_msg = || std::io::Error::other(format!("Invalid optional field {field}"));

        let (tag_text, rest) = field.split_once(':').ok_or_else(inv_opt_err_msg)?;
        let (type_text, string_value) = rest.split_once(':').ok_or_else(inv_opt_err_msg)?;

        let tag = SamOptField::parse_tag(tag_text)?;
        let type_code = SamOptField::parse_type(type_text)?;
        let opt_field = SamOptField::parse_value(tag, type_code, string_value)
            .with_context(|| format!("Failed to parse field '{field}'"))?;

        Ok(opt_field)
    }
}

impl AsRef<SamOptField> for SamOptField {
    fn as_ref(&self) -> &SamOptField {
        self
    }
}
