use crate::{data::err::ResultWithErrorContext, iter_utils::ProcessResultsExt, math::AnyInt, prelude::*};

/// Any optional fields stored in a SAM record, lazily parsed on an as-needed
/// basis.
///
/// Each optional field consists of a tag, value type, and value.
#[derive(Clone, Debug, Default)]
pub struct SamOptRaw(pub(super) Vec<String>);

impl SamOptRaw {
    /// Returns an empty collection of optional fields.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        SamOptRaw(Vec::new())
    }

    /// Returns [`SamOptRaw`] containing just a single field with the alignment
    /// score.
    ///
    /// The score is represented as an integer using the `AS` tag.
    #[inline]
    #[must_use]
    pub fn new_with_score<T: AnyInt + Into<i64>>(score: T) -> Self {
        let mut inner = Vec::with_capacity(1);
        inner.push(format!("AS:i:{score}", score = score.into()));
        SamOptRaw(inner)
    }

    /// Returns whether the optional data is empty.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of optional fields present.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Provides an iterator over the optional fields present (the tag names and
    /// parsed values).
    ///
    /// ## Limitations
    ///
    /// This iterator parses the fields lazily. If the [`SamOptRaw`] struct will
    /// be iterated over many times, consider parsing the fields once and
    /// collecting them.
    ///
    /// ## Errors
    ///
    /// The field in the SAM record must be of the form `TAG:TYPE:VALUE`. `TAG`
    /// cannot contain a colon. `TYPE` must be either `A`, `i`, `f`, `Z`, `H`,
    /// or `B`. `VALUE` must successfully parse into the corresponding type.
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = std::io::Result<SamOptField>> {
        self.as_view().iter()
    }

    /// Provides an iterator over the raw, unparsed optional fields present.
    #[inline]
    pub fn iter_raw(&self) -> std::slice::Iter<'_, String> {
        self.0.iter()
    }

    /// Returns the optional data for the provided tag, if it is present.
    ///
    /// ## Limitations
    ///
    /// This struct parses fields lazily and will repeat the computations each
    /// time [`get`] is called. The function runs in $O(n)$ time where $n$ is
    /// the number of fields. Validation of the field format is only performed
    /// where necessary.
    ///
    /// ## Errors
    ///
    /// The field in the SAM record must be of the form `TAG:TYPE:VALUE`. `TAG`
    /// cannot contain a colon. `TYPE` must be either `A`, `i`, `f`, `Z`, `H`,
    /// `B`. `VALUE` must successfully parse into the corresponding type.
    ///
    /// [`get`]: SamOptRaw::get
    pub fn get(&self, tag: &str) -> std::io::Result<Option<SamOptField>> {
        self.as_view().get(tag)
    }

    /// Adds an optional field to the [`SamOptRaw`] struct.
    ///
    /// ## Validity
    ///
    /// The tag name being pushed should not already be present in `self`.
    #[inline]
    pub fn push(&mut self, tag: &str, data: &SamOptValue) {
        self.0.push(format!("{tag}:{data}"));
    }
}

/// Any optional fields stored in a SAM record, lazily parsed on an as-needed
/// basis.
///
/// Each optional field consists of a tag, value type, and value.
#[derive(Copy, Clone, Debug, Default)]
pub struct SamOptRawView<'a>(pub(super) &'a [String]);

impl SamOptRawView<'_> {
    /// Returns an empty collection of optional fields.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        SamOptRawView(&[])
    }

    /// Returns whether the optional data is empty.
    #[inline]
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of optional fields present.
    #[inline]
    #[must_use]
    pub fn len(self) -> usize {
        self.0.len()
    }

    /// Provides an iterator over the optional fields present (the tag names and
    /// parsed values).
    ///
    /// ## Limitations
    ///
    /// This iterator parses the fields lazily. If the [`SamOptRaw`] struct will
    /// be iterated over many times, consider parsing the fields once and
    /// collecting them.
    ///
    /// ## Errors
    ///
    /// The field in the SAM record must be of the form `TAG:TYPE:VALUE`. `TAG`
    /// cannot contain a colon. `TYPE` must be either `A`, `i`, `f`, `Z`, `H`,
    /// or `B`. `VALUE` must successfully parse into the corresponding type.
    #[inline]
    pub fn iter(self) -> impl Iterator<Item = std::io::Result<SamOptField>> {
        self.0.iter().map(|field| {
            let inv_opt_err_msg = || std::io::Error::other(format!("Invalid optional field {field}"));

            let (tag_text, rest) = field.split_once(':').ok_or_else(inv_opt_err_msg)?;
            let (type_text, string_value) = rest.split_once(':').ok_or_else(inv_opt_err_msg)?;

            let tag = SamOptField::parse_tag(tag_text)?;
            let type_code = SamOptField::parse_type(type_text)?;
            let opt_field = SamOptField::parse_value(tag, type_code, string_value)
                .with_context(format!("Failed to parse field '{field}'"))?;
            Ok(opt_field)
        })
    }

    /// Provides an iterator over the raw, unparsed optional fields present.
    #[inline]
    pub fn iter_raw(&self) -> std::slice::Iter<'_, String> {
        self.0.iter()
    }

    /// Returns the optional data for the provided tag, if it is present.
    ///
    /// ## Limitations
    ///
    /// This struct parses fields lazily and will repeat the computations each
    /// time [`get`] is called. The function runs in $O(n)$ time where $n$ is
    /// the number of fields. Validation of the field format is only performed
    /// where necessary.
    ///
    /// ## Errors
    ///
    /// The field in the SAM record must be of the form `TAG:TYPE:VALUE`. `TAG`
    /// cannot contain a colon. `TYPE` must be either `A`, `i`, `f`, `Z`, `H`,
    /// `B`. `VALUE` must successfully parse into the corresponding type.
    ///
    /// [`get`]: SamOptRaw::get
    pub fn get(self, tag: &str) -> std::io::Result<Option<SamOptField>> {
        for field in self.0 {
            let inv_opt_err_msg = || std::io::Error::other(format!("Invalid optional field {field}"));
            let (this_tag, rest) = field.split_once(':').ok_or_else(inv_opt_err_msg)?;

            if this_tag == tag {
                let (type_text, string_value) = rest.split_once(':').ok_or_else(inv_opt_err_msg)?;

                let tag = SamOptField::parse_tag(this_tag)?;
                let type_code = SamOptField::parse_type(type_text)?;

                let opt_field = match SamOptField::parse_value(tag, type_code, string_value) {
                    Ok(opt_field) => opt_field,
                    Err(e) => {
                        return Err(std::io::Error::other(format!(
                            "Failed to parse field '{field}' due to error: {e}"
                        )));
                    }
                };
                return Ok(Some(opt_field));
            }
        }
        Ok(None)
    }
}

impl FromIterator<String> for SamOptRaw {
    /// Collects an iterator of strings into a [`SamOptRaw`] collection (each
    /// following the SAM file format for an optional field).
    ///
    /// ## Validity
    ///
    /// Each string should conform to the SAM file format for an optional field.
    /// Specifically, each field should be of the form `TAG:TYPE:VALUE`. `TAG`
    /// cannot contain a colon. `TYPE` must be either `A`, `i`, `f`, `Z`, `H`,
    /// or `B`. `VALUE` must successfully parse into the corresponding type.
    /// Furthermore, the tags should be unique.
    #[inline]
    fn from_iter<T: IntoIterator<Item = String>>(iter: T) -> Self {
        SamOptRaw(Vec::from_iter(iter))
    }
}

/// A parsed optional field in the SAM file format.
#[derive(Clone, Debug)]
pub struct SamOptField {
    /// The tag of the optional SAM field.
    pub tag:   [u8; 2],
    /// The value of the optional SAM field.
    pub value: SamOptValue,
}

impl SamOptField {
    /// Parses the tag for the optional SAM field from a string slice.
    fn parse_tag(tag: &str) -> std::io::Result<[u8; 2]> {
        let bytes = tag.as_bytes();
        if bytes.len() != 2 || !bytes[0].is_ascii_alphabetic() || !bytes[1].is_ascii_alphanumeric() {
            return Err(std::io::Error::other(format!("Invalid SAM optional tag: {tag}")));
        }
        Ok([bytes[0], bytes[1]])
    }

    /// Parses the type for the optional SAM field from a string slice.
    fn parse_type(type_text: &str) -> std::io::Result<char> {
        let mut type_chars = type_text.chars();
        let Some(typ) = type_chars.next() else {
            return Err(std::io::Error::other("Missing optional field type"));
        };
        if type_chars.next().is_some() {
            return Err(std::io::Error::other(format!("Invalid optional field type {type_text}")));
        }

        Ok(typ)
    }

    /// Parses a [`SamOptField`] from a tag, type, and value (as a string).
    ///
    /// ## Errors
    ///
    /// `type_code` must contain a valid character (`A`, `i`, `f`, `Z`, `H`, or
    /// `B`). The `string_value` must successfully parse into the corresponding
    /// type.
    fn parse_value(tag: [u8; 2], type_code: char, string_value: &str) -> std::io::Result<SamOptField> {
        match type_code {
            'A' => {
                let mut chars = string_value.chars();
                let Some(c) = chars.next() else {
                    return Err(std::io::Error::other("'A' field has empty value"));
                };
                if chars.next().is_some() {
                    return Err(std::io::Error::other("'A' field must contain exactly one character"));
                }
                if !c.is_ascii_graphic() {
                    return Err(std::io::Error::other("'A' field must be a printable ASCII character"));
                }
                Ok(SamOptField {
                    tag,
                    value: SamOptValue::Char(c as u8),
                })
            }
            'i' => {
                let parsed = string_value.parse::<i64>().with_context("Error parsing 'i' field")?;
                Ok(SamOptField {
                    tag,
                    value: SamOptValue::Int(parsed),
                })
            }
            'f' => {
                let parsed = string_value.parse::<f32>().with_context("Error parsing 'f' field")?;
                if !parsed.is_finite() {
                    return Err(std::io::Error::other("'f' field must be finite"));
                }
                Ok(SamOptField {
                    tag,
                    value: SamOptValue::Float(parsed),
                })
            }
            'Z' => {
                if !string_value.chars().all(|c| c == ' ' || c.is_ascii_graphic()) {
                    return Err(std::io::Error::other("'Z' field must contain printable ASCII characters"));
                }
                Ok(SamOptField {
                    tag,
                    value: SamOptValue::String(String::from(string_value)),
                })
            }
            'H' => {
                if !string_value.len().is_multiple_of(2) {
                    return Err(std::io::Error::other(format!(
                        "'H' field must contain an even number digits. Found {}",
                        string_value.len()
                    )));
                }
                if !string_value.as_bytes().iter().all(u8::is_ascii_hexdigit) {
                    return Err(std::io::Error::other("'H' field must contain hexadecimal digits"));
                }
                Ok(SamOptField {
                    tag,
                    value: SamOptValue::Hex(string_value.to_ascii_uppercase()),
                })
            }
            'B' => Ok(SamOptField {
                tag,
                value: SamOptValue::Array(
                    OptArray::parse_subtype_and_vals(string_value).with_context("Failed to parse 'B' array")?,
                ),
            }),
            _ => Err(std::io::Error::other(format!(
                "Unsupported SAM optional field type {type_code}"
            ))),
        }
    }

    /// Returns the stored character from the [`SamOptField`], or [`None`] if a
    /// different variant is present.
    #[inline]
    #[must_use]
    pub fn char(self) -> Option<u8> {
        match self.value {
            SamOptValue::Char(c) => Some(c),
            _ => None,
        }
    }

    /// Returns the stored integer from the [`SamOptField`], or [`None`] if a
    /// different variant is present.
    #[inline]
    #[must_use]
    pub fn int(self) -> Option<i64> {
        match self.value {
            SamOptValue::Int(i) => Some(i),
            _ => None,
        }
    }

    /// Returns the stored floating point number from the [`SamOptField`], or
    /// [`None`] if a different variant is present.
    #[inline]
    #[must_use]
    pub fn float(self) -> Option<f32> {
        match self.value {
            SamOptValue::Float(f) => Some(f),
            _ => None,
        }
    }

    /// Returns the stored string from the [`SamOptField`], or [`None`] if a
    /// different variant is present.
    #[inline]
    #[must_use]
    pub fn string(self) -> Option<String> {
        match self.value {
            SamOptValue::String(f) => Some(f),
            _ => None,
        }
    }

    /// Returns the stored hex string from the [`SamOptField`], or [`None`] if a
    /// different variant is present.
    ///
    /// For example, the six-character hex string "1AE301" represents the byte
    /// array `[0x1a, 0xe3, 0x01]`.
    #[inline]
    #[must_use]
    pub fn hex(self) -> Option<String> {
        match self.value {
            SamOptValue::Hex(f) => Some(f),
            _ => None,
        }
    }

    /// Returns the stored [`OptArray`] from the [`SamOptField`], or [`None`] if
    /// a different variant is present.
    #[inline]
    #[must_use]
    pub fn array(self) -> Option<OptArray> {
        match self.value {
            SamOptValue::Array(f) => Some(f),
            _ => None,
        }
    }
}

/// The value of an optional field (for the SAM file format).
#[derive(Clone, PartialEq, Debug)]
pub enum SamOptValue {
    /// A printable character (type code `A`).
    Char(u8),
    /// A signed integer (type code `i`).
    Int(i64),
    /// A single-precision floating number (type code `f`).
    Float(f32),
    /// A printable string, including space (type code `Z`).
    String(String),
    /// A byte array in the hex format (type code `H`).
    ///
    /// For example, the six-character hex string "1AE301" represents the byte
    /// array `[0x1a, 0xe3, 0x01]`.
    Hex(String),
    // An integer or numeric array (type code `B`).
    Array(OptArray),
}

/// The data array for a [`SamOptValue::Array`] variant (type code `B`).
#[derive(Clone, PartialEq, Debug)]
pub enum OptArray {
    /// Array subtype code `c`.
    I8(Vec<i8>),
    /// Array subtype code `C`.
    U8(Vec<u8>),
    /// Array subtype code `s`.
    I16(Vec<i16>),
    /// Array subtype code `S`.
    U16(Vec<u16>),
    /// Array subtype code `i`.
    I32(Vec<i32>),
    /// Array subtype code `I`.
    U32(Vec<u32>),
    /// Array subtype code `f`.
    F32(Vec<f32>),
}

impl OptArray {
    /// Parses the array of optional SAM fields with type `B`.
    ///
    /// ## Errors
    ///
    /// The first letter in the array indicates the type of numbers in the
    /// following comma-separated array. The letter can be one of `c`, `C`, `s`,
    /// `S`, `i`, `I`, or `f`.
    fn parse_subtype_and_vals(string_value: &str) -> std::io::Result<Self> {
        let mut pieces = string_value.split(',');
        let Some(subtype) = pieces.next() else {
            return Err(std::io::Error::other("Missing subtype"));
        };

        match subtype {
            "c" => {
                let values = pieces
                    .map(str::parse::<i8>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'c' subtype (`i8`)")?;

                Ok(OptArray::I8(values))
            }
            "C" => {
                let values = pieces
                    .map(str::parse::<u8>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'C' subtype (`u8`)")?;
                Ok(OptArray::U8(values))
            }
            "s" => {
                let values = pieces
                    .map(str::parse::<i16>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 's' subtype (`i16`)")?;
                Ok(OptArray::I16(values))
            }
            "S" => {
                let values = pieces
                    .map(str::parse::<u16>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'S' subtype (`u16`)")?;
                Ok(OptArray::U16(values))
            }
            "i" => {
                let values = pieces
                    .map(str::parse::<i32>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'i' subtype (`i32`)")?;
                Ok(OptArray::I32(values))
            }
            "I" => {
                let values = pieces
                    .map(str::parse::<u32>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'I' subtype (`u32`)")?;
                Ok(OptArray::U32(values))
            }
            "f" => {
                let values = pieces
                    .map(|value| {
                        let parsed = value.parse::<f32>().map_err(std::io::Error::other)?;
                        if !parsed.is_finite() {
                            return Err(std::io::Error::other("'B:f' values must be finite"));
                        }
                        Ok(parsed)
                    })
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'f' subtype (`f32`)")?;
                Ok(OptArray::F32(values))
            }
            _ => Err(std::io::Error::other(format!("Unsupported subtype {subtype}"))),
        }
    }
}
