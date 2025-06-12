use crate::{
    alignment::phmm::{
        CorePhmm, EmissionParams, GlobalPhmm, LayerParams, PhmmNumber,
        PhmmState::{self, *},
        TransitionParams,
    },
    data::{
        ByteIndexMap,
        mappings::{AA_UNAMBIG_PROFILE_MAP, DNA_UNAMBIG_PROFILE_MAP},
    },
    unwrap_or_return_some_err,
};
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Error as IOError, ErrorKind, Lines, Write},
    marker::PhantomData,
    path::Path,
};

/// A struct providing methods for parsing a [`GlobalPhmm`] from a SAM (sequence
/// alignment and modeling system) model file, using `f32` to represent all
/// model parameters.
pub struct SamHmmParser;

impl SamHmmParser {
    /// Parse a DNA pHMM from a SAM model file.
    ///
    /// ## Errors
    ///
    /// * The file must meet the specifications for a SAM model file
    /// * The file must be a model file (not a regularizer or null model)
    /// * The specified alphabet must be dna
    /// * No negative indices are allowed in layer names
    pub fn parse_dna_model(filename: impl AsRef<Path>) -> std::io::Result<GlobalPhmm<f32, 4>> {
        SupportedConfig::parse_sam_model_file(filename)
    }

    /// Parse a protein pHMM from a SAM model file.
    ///
    /// ## Errors
    ///
    /// * The file must meet the specifications for a SAM model file
    /// * The file must be a model file (not a regularizer or null model)
    /// * The specified alphabet must be protein
    /// * No negative indices are allowed in layer names
    pub fn parse_protein_model(filename: impl AsRef<Path>) -> Result<GlobalPhmm<f32, 20>, std::io::Error> {
        SupportedConfig::parse_sam_model_file(filename)
    }
}

/// A struct providing methods for parsing a [`GlobalPhmm`] from a SAM (sequence
/// alignment and modeling system) model file, using type `T` to represent all
/// model parameters.
struct GenericSamHmmParser<T>(PhantomData<T>);

impl<T: PhmmNumber> GenericSamHmmParser<T> {
    /// Parse a DNA pHMM from a SAM model file.
    ///
    /// ## Errors
    ///
    /// * The file must meet the specifications for a SAM model file
    /// * The file must be a model file (not a regularizer or null model)
    /// * The specified alphabet must be dna
    /// * No negative indices are allowed in layer names
    fn parse_dna_model(filename: impl AsRef<Path>) -> Result<GlobalPhmm<T, 4>, std::io::Error> {
        SupportedConfig::parse_sam_model_file(filename)
    }

    /// Parse a protein pHMM from a SAM model file.
    ///
    /// ## Errors
    ///
    /// * The file must meet the specifications for a SAM model file
    /// * The file must be a model file (not a regularizer or null model)
    /// * The specified alphabet must be protein
    /// * No negative indices are allowed in layer names
    fn parse_protein_model(filename: impl AsRef<Path>) -> Result<GlobalPhmm<T, 20>, std::io::Error> {
        SupportedConfig::parse_sam_model_file(filename)
    }
}

/// A struct providing methods for writing a [`GlobalPhmm`] to a SAM (sequence
/// alignment and modeling system) model file.
pub struct SamHmmWriter;

impl SamHmmWriter {
    /// Write a DNA pHMM to a SAM model file.
    ///
    /// ## Errors
    ///
    /// * IO errors (when creating file or writing to it)
    /// * The mapping of the pHMM must correspond to DNA
    /// * The model must have at least one layer
    #[inline]
    pub fn write_dna_model<T: PhmmNumber>(filename: impl AsRef<Path>, model: &GlobalPhmm<T, 4>) -> std::io::Result<()> {
        SupportedConfig::write_sam_model_file(filename, model)
    }

    /// Write a protein pHMM to a SAM model file.
    ///
    /// ## Errors
    ///
    /// * IO errors (when creating file or writing to it)
    /// * The mapping of the pHMM must correspond to DNA
    /// * The model must have at least one layer
    #[inline]
    pub fn write_protein_model<T: PhmmNumber>(filename: impl AsRef<Path>, model: &GlobalPhmm<T, 20>) -> std::io::Result<()> {
        SupportedConfig::write_sam_model_file(filename, model)
    }
}

/// A list of the pHMM transition probabilities to copy to the previous layer
/// when parsing.
const PARAMS_TO_COPY_PARSING: [(PhmmState, PhmmState); 6] = [
    (Delete, Delete),
    (Delete, Match),
    (Match, Delete),
    (Match, Match),
    (Insert, Delete),
    (Insert, Match),
];

/// A zero-size struct through which supported parser configurations can be
/// specified (via implementing the [`SamHmmConfig`] trait).
struct SupportedConfig;

/// A trait providing the ability to handle different alphabets in SAM files.
trait SamHmmConfig<const S: usize, const L: usize> {
    /// Parses the alphabet line to the appropriate [`ByteIndexMap`]. The input
    /// is the content in the alphabet line of the model file, converted to
    /// uppercase with the "alphabet" prefix and any leading/trailing whitespace
    /// removed.
    fn parse_mapping(mapping: &str) -> std::io::Result<&'static ByteIndexMap<S>>;

    /// Converts a [`ByteIndexMap`] to the appropriate alphabet name for the SAM
    /// file.
    fn unparse_mapping(mapping: &'static ByteIndexMap<S>) -> std::io::Result<&'static str>;

    /// Given a flat array of `L=2*S+9` parameters, converts them into
    /// [`LayerParams`].
    fn group_params<T: PhmmNumber>(params: [T; L]) -> LayerParams<T, S>;

    /// The inverse of `group_params`, taking a layer and flattening it to SAM
    /// parameters.
    fn ungroup_params<T: PhmmNumber>(params: &LayerParams<T, S>) -> [T; L];

    /// Parse a SAM model file representing a pHMM.
    ///
    /// ## Errors
    ///
    /// * The file must meet the specifications for a SAM model file
    /// * The file must be a model file (not a regularizer or null model)
    /// * The specified alphabet must be compatible with the selected `S` and
    ///   `L`
    /// * No negative indices are allowed in layer names
    /// * At least three layers must be present in the model
    fn parse_sam_model_file<T: PhmmNumber, P>(filename: P) -> Result<GlobalPhmm<T, S>, IOError>
    where
        P: AsRef<Path>,
        SupportedConfig: SamHmmConfig<S, L>, {
        let mut lines = LineIterator::new(filename)?;

        // Read initial line
        validate_model_line(&mut lines)?;

        // Extract and read alphabet line
        let mapping = Self::parse_alphabet_line(&mut lines)?;

        // Read all model layers
        let mut layers = LayerIter::<T, S, L>::new(&mut lines)?.collect::<Result<Vec<_>, IOError>>()?;

        let [first_layer, .., last_layer] = layers.as_mut_slice() else {
            return Err(IOError::new(
                ErrorKind::InvalidData,
                "At least three layers must be specified in the SAM model file",
            ));
        };

        first_layer.transition[(Delete, Delete)] = T::INFINITY;
        first_layer.transition[(Delete, Match)] = T::INFINITY;
        first_layer.transition[(Delete, Insert)] = T::INFINITY;

        last_layer.transition[(Delete, Delete)] = T::INFINITY;
        last_layer.transition[(Insert, Delete)] = T::INFINITY;
        last_layer.transition[(Match, Delete)] = T::INFINITY;

        Ok(GlobalPhmm {
            mapping,
            // Validity: We already verified `layers` has at least two entries
            core: CorePhmm::new_unchecked(layers),
        })
    }

    /// Write a global pHMM to a file, following the SAM format.
    ///
    /// No comment lines are included.
    ///
    /// ## Errors
    ///
    /// * IO errors (when creating file or writing to it)
    /// * The mapping of the pHMM must correspond to DNA
    /// * The model must have at least one layer
    fn write_sam_model_file<T: PhmmNumber, P>(filename: P, model: &GlobalPhmm<T, S>) -> std::io::Result<()>
    where
        P: AsRef<Path>,
        SupportedConfig: SamHmmConfig<S, L>, {
        if model.core.0.len() < 2 {
            return Err(IOError::new(
                ErrorKind::InvalidData,
                "At least two layers must be present in the model!",
            ));
        }

        let mut writer = BufWriter::new(File::create(filename.as_ref())?);

        writeln!(writer, "MODEL")?;
        writeln!(writer, "alphabet {}", Self::unparse_mapping(model.mapping)?)?;

        let [first_layer, rest @ ..] = model.core.0.as_slice() else {
            return Err(IOError::new(
                ErrorKind::InvalidData,
                "At least two layers must be present in the model!",
            ));
        };

        let mut current_layer = LayerParams::<T, S>::default();

        write!(writer, "0 ")?;
        current_layer.transition[Insert] = first_layer.transition[Insert];
        current_layer.emission_insert = first_layer.emission_insert.clone();
        print_params(&mut writer, Self::ungroup_params(&current_layer))?;
        current_layer = first_layer.clone();

        let mut i = 1;
        for layer in rest {
            write!(writer, "{i} ")?;
            current_layer.transition[Insert] = layer.transition[Insert];
            current_layer.emission_insert = layer.emission_insert.clone();
            print_params(&mut writer, Self::ungroup_params(&current_layer))?;
            current_layer = layer.clone();
            i += 1;
        }

        write!(writer, "END ")?;
        print_params(&mut writer, Self::ungroup_params(&current_layer))?;
        writeln!(writer)?;
        writeln!(writer, "ENDMODEL")
    }

    /// Parse the alphabet line of the model. This consumes one line of the file
    fn parse_alphabet_line(lines: &mut LineIterator) -> Result<&'static ByteIndexMap<S>, IOError> {
        let line = lines
            .next()
            .ok_or(IOError::new(ErrorKind::InvalidData, "Could not locate alphabet line in file"))??
            .to_ascii_uppercase();

        let mut tokens = line.split_whitespace();

        let field = tokens.next().unwrap_or("");
        if !field.eq_ignore_ascii_case("ALPHABET") {
            return Err(IOError::new(ErrorKind::InvalidData, "Expected alphabet line"));
        }

        let Some(mut alphabet) = tokens.next() else {
            return Err(IOError::new(ErrorKind::InvalidData, "Alphabet was missing"));
        };
        alphabet = alphabet.trim();

        Self::parse_mapping(alphabet)
    }

    /// Parses the parameters for a pHMM layer from a set of tokens
    /// (`rest_of_line`) followed by as many lines from `lines` as needed. There
    /// are assumed to be `L=2*S+9` parameters, representing the 9 transition
    /// probabilities, `S` match emission probabilities, and `S` insert emission
    /// probabilities.
    ///
    /// ## Errors
    ///
    /// * `rest_of_line` and `lines` must contain sufficiently many tokens to
    ///   fill a full set of parameters
    /// * All parameters must parse successfully, and there should be no extra
    ///   parameters on any line (see [`fill_params_from_iter`])
    ///
    /// [`fill_params_from_iter`]: crate::alignment::phmm::sam_parser::SamHmmConfig::fill_params_from_iter
    fn parse_layer_params<'a, T: PhmmNumber>(
        mut rest_of_line: impl Iterator<Item = &'a str>, lines: &mut LineIterator,
    ) -> Result<LayerParams<T, S>, IOError> {
        let mut params = [T::ZERO; L];

        let mut i = Self::fill_params_from_iter::<T>(rest_of_line.by_ref(), &mut params, 0)?;

        if i == L {
            return Ok(Self::group_params(params));
        }

        for line in lines {
            let line = line?;
            let mut split = line.split_whitespace();
            i = Self::fill_params_from_iter::<T>(split.by_ref(), &mut params, i)?;
            if i == L {
                break;
            }
        }

        if i < params.len() {
            return Err(IOError::new(
                ErrorKind::InvalidData,
                "The file ended before a full set of parameters was parsed",
            ));
        }

        Ok(Self::group_params(params))
    }

    /// Uses the parameters in `iter` to fill spots in `params`, starting with
    /// `i` and continuing until the end of `params` or `iter`. Returns the
    /// index of the next spot to fill (or `params.len()` if it is full).
    ///
    /// ## Errors
    ///
    /// * The parameters must successfully parse (see [`parse_param`])
    /// * If `params` gets filled, then `iter` must not contain any more
    ///   elements
    fn fill_params_from_iter<'a, T: PhmmNumber>(
        iter: &mut impl Iterator<Item = &'a str>, params: &mut [T], mut i: usize,
    ) -> Result<usize, std::io::Error> {
        for param in iter.take(L - i) {
            params[i] = parse_param(param)?;
            i += 1;
        }
        if iter.next().is_some() {
            return Err(IOError::new(
                ErrorKind::InvalidData,
                "Layers of the model must be separated by line breaks",
            ));
        }
        Ok(i)
    }
}

impl SamHmmConfig<4, 17> for SupportedConfig {
    #[inline]
    fn parse_mapping(mapping: &str) -> Result<&'static ByteIndexMap<4>, IOError> {
        match mapping {
            "DNA" => Ok(&DNA_UNAMBIG_PROFILE_MAP),
            "PROTEIN" => Err(IOError::new(
                ErrorKind::InvalidData,
                "A protein alphabet was found, but a DNA alphabet was expected",
            )),
            _ => Err(IOError::new(ErrorKind::InvalidData, "Unsupported alphabet specified")),
        }
    }

    #[inline]
    fn unparse_mapping(mapping: &'static ByteIndexMap<4>) -> std::io::Result<&'static str> {
        if mapping == &DNA_UNAMBIG_PROFILE_MAP {
            return Ok("DNA");
        }
        Err(IOError::new(ErrorKind::InvalidData, "The mapping is unsupported by SAM!"))
    }

    #[inline]
    fn group_params<T: PhmmNumber>(params: [T; 17]) -> LayerParams<T, 4> {
        let transition = TransitionParams([
            [params[0], params[1], params[2]],
            [params[3], params[4], params[5]],
            [params[6], params[7], params[8]],
        ]);

        let emission_match = EmissionParams([params[9], params[11], params[10], params[12]]);
        let emission_insert = EmissionParams([params[13], params[15], params[14], params[16]]);

        LayerParams {
            transition,
            emission_match,
            emission_insert,
        }
    }

    #[inline]
    fn ungroup_params<T: PhmmNumber>(params: &LayerParams<T, 4>) -> [T; 17] {
        let mut out = [T::ZERO; 17];
        out[0..3].copy_from_slice(&params.transition.0[0]);
        out[3..6].copy_from_slice(&params.transition.0[1]);
        out[6..9].copy_from_slice(&params.transition.0[2]);
        out[9] = params.emission_match[0];
        out[10] = params.emission_match[2];
        out[11] = params.emission_match[1];
        out[12] = params.emission_match[3];
        out[13] = params.emission_insert[0];
        out[14] = params.emission_insert[2];
        out[15] = params.emission_insert[1];
        out[16] = params.emission_insert[3];
        out
    }
}

impl SamHmmConfig<20, 49> for SupportedConfig {
    #[inline]
    fn parse_mapping(mapping: &str) -> Result<&'static ByteIndexMap<20>, IOError> {
        match mapping {
            "DNA" => Err(IOError::new(
                ErrorKind::InvalidData,
                "A DNA alphabet was found, but a protein alphabet was expected",
            )),
            "PROTEIN" => Ok(&AA_UNAMBIG_PROFILE_MAP),
            _ => Err(IOError::new(ErrorKind::InvalidData, "Unsupported alphabet specified")),
        }
    }

    #[inline]
    fn unparse_mapping(mapping: &'static ByteIndexMap<20>) -> std::io::Result<&'static str> {
        if mapping == &AA_UNAMBIG_PROFILE_MAP {
            return Ok("PROTEIN");
        }
        Err(IOError::new(ErrorKind::InvalidData, "The mapping is unsupported by SAM!"))
    }

    #[inline]
    fn group_params<T: PhmmNumber>(params: [T; 49]) -> LayerParams<T, 20> {
        let (transition, rest) = params.split_at(9);

        let transition = TransitionParams([
            [transition[0], transition[1], transition[2]],
            [transition[3], transition[4], transition[5]],
            [transition[6], transition[7], transition[8]],
        ]);
        let (emission_match, emission_insert) = rest.split_at(20);

        let emission_match = EmissionParams(emission_match.try_into().unwrap());
        let emission_insert = EmissionParams(emission_insert.try_into().unwrap());

        LayerParams {
            transition,
            emission_match,
            emission_insert,
        }
    }

    #[inline]
    fn ungroup_params<T: PhmmNumber>(params: &LayerParams<T, 20>) -> [T; 49] {
        let mut out = [T::ZERO; 49];
        out[0..3].copy_from_slice(&params.transition.0[0]);
        out[3..6].copy_from_slice(&params.transition.0[1]);
        out[6..9].copy_from_slice(&params.transition.0[2]);
        out[9..29].copy_from_slice(&params.emission_match.0);
        out[29..49].copy_from_slice(&params.emission_insert.0);
        out
    }
}

/// Parse a layer name. BEGIN returns 0, END returns None, otherwise it is
/// parsed as a usize and returned.
///
/// ## Errors
///
/// * Negatively-numbered nodes are not supported
/// * Any layer name other than BEGIN or END must successfully parse to a usize
fn parse_layer_name(token: &str) -> Result<Option<usize>, IOError> {
    if token.eq_ignore_ascii_case("BEGIN") {
        Ok(Some(0))
    } else if token.eq_ignore_ascii_case("END") {
        Ok(None)
    } else if token.starts_with('-') {
        Err(IOError::new(
            ErrorKind::InvalidData,
            "Negatively-numbered nodes in models are not supported. Consider using a prior version of SAM's hmmconvert to convert the model",
        ))
    } else if let Ok(layer) = token.parse::<usize>() {
        Ok(Some(layer))
    } else {
        Err(IOError::new(ErrorKind::InvalidData, "Could not parse the layer name"))
    }
}

/// Extracts and validates the "MODEL" line of the file
#[inline]
fn validate_model_line(lines: &mut LineIterator) -> Result<(), IOError> {
    let line = lines
        .next()
        .ok_or(IOError::new(ErrorKind::InvalidData, "Could not locate initial line in file"))??;

    let token = line.split_whitespace().next().unwrap_or("");

    match token {
        "MODEL" => Ok(()),
        "REGULARIZER" => Err(IOError::new(
            ErrorKind::InvalidData,
            "REGULARIZER is not supported by this parser",
        )),
        "NULLMODEL" => Err(IOError::new(
            ErrorKind::InvalidData,
            "NULLMODEL is not supported by this parser",
        )),
        _ => Err(IOError::new(ErrorKind::InvalidData, "Could not locate initial line in file")),
    }
}

/// Parse a parameter from a string, and convert it to its negative natural
/// logarithm. A parameter that is 0 is converted to $\infty$.
///
/// ## Errors
///
/// * Negative parameters are not allowed
/// * The parameter must succeed when parsing to type `T`
#[inline]
fn parse_param<T: PhmmNumber>(prob: &str) -> Result<T, IOError> {
    if prob.starts_with('-') {
        return Err(IOError::new(
            ErrorKind::InvalidData,
            format!("A negative parameter was found: {prob}"),
        ));
    }
    let prob = prob.parse::<f64>().map_err(|_| {
        IOError::new(
            ErrorKind::InvalidData,
            format!("When parsing the model parameters, the value {prob} could not be parsed as a float"),
        )
    })?;
    Ok(T::from_prob(prob))
}

/// An iterator over the layers in a SAM pHMM file, rearranged to be compatible
/// with Zoe.
struct LayerIter<'a, T, const S: usize, const L: usize> {
    /// The underlying raw layers of the SAM pHMM file
    raw_layers: RawLayerIter<'a, T, S, L>,
    /// The last layer that was parsed, but not yet yielded
    last_layer: LayerParams<T, S>,
}

impl<'a, T: PhmmNumber, const S: usize, const L: usize> LayerIter<'a, T, S, L>
where
    SupportedConfig: SamHmmConfig<S, L>,
{
    /// Creates a new iterator over the layers in a SAM pHMM file, rearranged to
    /// be compatible with Zoe.
    ///
    /// ## Errors
    ///
    /// * IO errors
    /// * At least one model layer must be present
    fn new(lines: &'a mut LineIterator) -> Result<Self, IOError> {
        let mut raw_layers = RawLayerIter::new(lines);
        let last_layer = raw_layers
            .next()
            .ok_or(IOError::new(ErrorKind::InvalidData, "No model layers found!"))??;
        Ok(Self { raw_layers, last_layer })
    }
}

impl<T: PhmmNumber, const S: usize, const L: usize> Iterator for LayerIter<'_, T, S, L>
where
    SupportedConfig: SamHmmConfig<S, L>,
{
    type Item = Result<LayerParams<T, S>, IOError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut next_layer = unwrap_or_return_some_err!(self.raw_layers.next()?);

        for param in PARAMS_TO_COPY_PARSING {
            self.last_layer.transition[param] = next_layer.transition[param];
        }
        self.last_layer.emission_match = std::mem::replace(&mut next_layer.emission_match, EmissionParams([T::ZERO; S]));

        Some(Ok(std::mem::replace(&mut self.last_layer, next_layer)))
    }
}

/// An iterator over the "raw" layers in a SAM pHMM file. These are the layers
/// as grouped in the model file, which is a different grouping than how
/// [`GlobalPhmm`] represents them.
struct RawLayerIter<'a, T, const S: usize, const L: usize> {
    /// The underlying iterator of lines in the file
    lines:          &'a mut LineIterator,
    /// The next layer number expected by the model. This is `None` when
    /// ENDMODEL has been parsed
    expected_layer: Option<usize>,
    phantom:        PhantomData<T>,
}

impl<'a, T, const S: usize, const L: usize> RawLayerIter<'a, T, S, L> {
    /// Creates a new [`RawLayerIter`] from an iterator of the lines in the
    /// file.
    fn new(lines: &'a mut LineIterator) -> Self {
        Self {
            lines,
            expected_layer: Some(0),
            phantom: PhantomData,
        }
    }
}

impl<T: PhmmNumber, const S: usize, const L: usize> Iterator for RawLayerIter<'_, T, S, L>
where
    SupportedConfig: SamHmmConfig<S, L>,
{
    type Item = Result<LayerParams<T, S>, IOError>;

    fn next(&mut self) -> Option<Self::Item> {
        // Get the next expected layer, which is `None` only if we have already
        // parsed ENDMODEL
        let expected_layer = self.expected_layer?;

        // Get line containing layer number/name
        let Some(line) = unwrap_or_return_some_err!(self.lines.next().transpose()) else {
            return Some(Err(IOError::new(ErrorKind::InvalidData, "File ended before ENDMODEL")));
        };

        // Get token corresponding to layer name/number
        let mut tokens = line.split_whitespace();
        let token = tokens.next().unwrap_or("");

        if token.eq_ignore_ascii_case("ENDMODEL") {
            self.expected_layer = None;
            return None;
        }

        if token.eq_ignore_ascii_case("FREQAVE") {
            unwrap_or_return_some_err!(SupportedConfig::parse_layer_params::<T>(tokens, self.lines));
            return self.next();
        }

        // Parse the layer name/number
        if let Some(layer_number) = unwrap_or_return_some_err!(parse_layer_name(token)) {
            if layer_number != expected_layer {
                return Some(Err(IOError::new(
                    ErrorKind::InvalidData,
                    format!("Found model node {layer_number}, expected model node {expected_layer}"),
                )));
            }

            self.expected_layer = Some(expected_layer + 1);
        } else {
            // Reached END, so we expect no more layers
            self.expected_layer = None;
        }

        // Parse the layer
        let layer = unwrap_or_return_some_err!(SupportedConfig::parse_layer_params(tokens, self.lines));

        // If this was the END layer, ensure we also parse ENDMODEL
        if self.expected_layer.is_none() {
            if let Some(line) = self.lines.next()
                && let Some(token) = unwrap_or_return_some_err!(line).split_whitespace().next()
            {
                if !token.eq_ignore_ascii_case("ENDMODEL") {
                    return Some(Err(IOError::new(
                        ErrorKind::InvalidData,
                        format!("Unexpected token {token} found between END layer and ENDMODEL"),
                    )));
                }
            } else {
                return Some(Err(IOError::new(ErrorKind::InvalidData, "Failed to find ENDMODEL")));
            }
        }

        Some(Ok(layer))
    }
}

/// An iterator over the lines in a SAM file, skipping empty lines, lines
/// containing only whitespace, and comment lines.
struct LineIterator {
    lines: Lines<BufReader<File>>,
}

impl LineIterator {
    /// Create a new [`LineIterator`] object from a file. A `BufReader` is
    /// automatically used.
    #[inline]
    fn new<P: AsRef<Path>>(filename: P) -> Result<Self, IOError> {
        Ok(Self {
            lines: BufReader::new(File::open(filename)?).lines(),
        })
    }
}

impl Iterator for LineIterator {
    type Item = Result<String, IOError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line = match self.lines.next()? {
                Ok(line) => line,
                Err(e) => return Some(Err(e)),
            };

            let line = line.trim().to_string();

            if !line.is_empty() && !line.starts_with('%') {
                return Some(Ok(line));
            }
        }
    }
}

/// Print the parameters for a layer to a single line, space-separated, and add
/// a newline.
#[inline]
fn print_params<T: PhmmNumber, const N: usize>(writer: &mut impl Write, params: [T; N]) -> std::io::Result<()> {
    let mut params = params.into_iter().map(super::PhmmNumber::to_prob::<f32>);
    let Some(param) = params.next() else { return Ok(()) };
    write!(writer, "{param:.6}")?;
    for param in params {
        write!(writer, " {param:.6}")?;
    }
    writeln!(writer)
}
