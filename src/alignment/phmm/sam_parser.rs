use crate::{
    alignment::phmm::{
        EmissionParams, LayerParams, Phmm,
        PhmmState::{self, *},
        TransitionParams,
    },
    data::{
        ByteIndexMap,
        mappings::{AA_UNAMBIG_PROFILE_MAP, DNA_UNAMBIG_PROFILE_MAP},
    },
    math::Float,
    unwrap_or_return_some_err,
};
use std::{
    fs::File,
    io::{BufRead, BufReader, Error as IOError, ErrorKind, Lines},
    marker::PhantomData,
    path::Path,
};

/// A struct providing methods for parsing a [`Phmm`] from a SAM (sequence
/// alignment and modeling system) model file, using `f32` to represent all
/// model parameters.
pub struct SamParser;

impl SamParser {
    /// Parse a DNA pHMM from a SAM model file.
    ///
    /// ## Errors
    ///
    /// * The file must meet the specifications for a SAM model file
    /// * The file must be a model file (not a regularizer or null model)
    /// * The specified alphabet must be dna
    /// * No negative indices are allowed in layer names
    pub fn parse_dna_model(filename: impl AsRef<Path>) -> Result<Phmm<f32, 4>, std::io::Error> {
        SupportedConfig::parse_sam_file(filename)
    }

    /// Parse a protein pHMM from a SAM model file.
    ///
    /// ## Errors
    ///
    /// * The file must meet the specifications for a SAM model file
    /// * The file must be a model file (not a regularizer or null model)
    /// * The specified alphabet must be protein
    /// * No negative indices are allowed in layer names
    pub fn parse_protein_model(filename: impl AsRef<Path>) -> Result<Phmm<f32, 20>, std::io::Error> {
        SupportedConfig::parse_sam_file(filename)
    }
}

/// A struct providing methods for parsing a [`Phmm`] from a SAM (sequence
/// alignment and modeling system) model file, using type `T` to represent all
/// model parameters.
struct GenericSamParser<T>(PhantomData<T>);

impl<T: Float> GenericSamParser<T> {
    /// Parse a DNA pHMM from a SAM model file.
    ///
    /// ## Errors
    ///
    /// * The file must meet the specifications for a SAM model file
    /// * The file must be a model file (not a regularizer or null model)
    /// * The specified alphabet must be dna
    /// * No negative indices are allowed in layer names
    fn parse_dna_model(filename: impl AsRef<Path>) -> Result<Phmm<T, 4>, std::io::Error> {
        SupportedConfig::parse_sam_file(filename)
    }

    /// Parse a protein pHMM from a SAM model file.
    ///
    /// ## Errors
    ///
    /// * The file must meet the specifications for a SAM model file
    /// * The file must be a model file (not a regularizer or null model)
    /// * The specified alphabet must be protein
    /// * No negative indices are allowed in layer names
    fn parse_protein_model(filename: impl AsRef<Path>) -> Result<Phmm<T, 20>, std::io::Error> {
        SupportedConfig::parse_sam_file(filename)
    }
}

/// A list of the pHMM transition probabilities to copy to the previous layer
/// when parsing.
const PARAMS_TO_COPY: [(PhmmState, PhmmState); 6] = [
    (Delete, Delete),
    (Delete, Match),
    (Match, Delete),
    (Match, Match),
    (Insert, Delete),
    (Insert, Match),
];

/// A zero-size struct through which supported parser configurations can be
/// specified (via implementing the [`ParserConfig`] trait).
struct SupportedConfig;

/// A trait providing the ability to parse different alphabets from SAM files.
trait ParserConfig<const S: usize, const L: usize> {
    /// Retrieves the appropriate [`ByteIndexMap`]. The input is the content in
    /// the alphabet line of the model file, converted to uppercase with the
    /// "alphabet" prefix and any leading/trailing whitespace removed.
    fn get_mapping(mapping: &str) -> Result<&'static ByteIndexMap<S>, IOError>;

    /// Given a flat array of `L=2*S+9` parameters, converts them into
    /// [`LayerParams`].
    fn group_params<T: Float>(params: [T; L]) -> LayerParams<T, S>;

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
    fn parse_sam_file<T: Float, P>(filename: P) -> Result<Phmm<T, S>, IOError>
    where
        P: AsRef<Path>,
        SupportedConfig: ParserConfig<S, L>, {
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
                "At least three layers must be specified",
            ));
        };

        first_layer.transition[(Delete, Delete)] = T::INFINITY;
        first_layer.transition[(Delete, Match)] = T::INFINITY;
        first_layer.transition[(Delete, Insert)] = T::INFINITY;

        last_layer.transition[(Delete, Delete)] = T::INFINITY;
        last_layer.transition[(Insert, Delete)] = T::INFINITY;
        last_layer.transition[(Match, Delete)] = T::INFINITY;

        Ok(Phmm { mapping, params: layers })
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

        Self::get_mapping(alphabet)
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
    fn parse_layer_params<'a, T: Float>(
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
    fn fill_params_from_iter<'a, T: Float>(
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

impl ParserConfig<4, 17> for SupportedConfig {
    fn get_mapping(mapping: &str) -> Result<&'static ByteIndexMap<4>, IOError> {
        match mapping {
            "DNA" => Ok(&DNA_UNAMBIG_PROFILE_MAP),
            "PROTEIN" => Err(IOError::new(
                ErrorKind::InvalidData,
                "A protein alphabet was found, but a DNA alphabet was expected",
            )),
            _ => Err(IOError::new(ErrorKind::InvalidData, "Unsupported alphabet specified")),
        }
    }

    fn group_params<T: Float>(params: [T; 17]) -> LayerParams<T, 4> {
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
}

impl ParserConfig<20, 49> for SupportedConfig {
    fn get_mapping(mapping: &str) -> Result<&'static ByteIndexMap<20>, IOError> {
        match mapping {
            "DNA" => Err(IOError::new(
                ErrorKind::InvalidData,
                "A DNA alphabet was found, but a protein alphabet was expected",
            )),
            "PROTEIN" => Ok(&AA_UNAMBIG_PROFILE_MAP),
            _ => Err(IOError::new(ErrorKind::InvalidData, "Unsupported alphabet specified")),
        }
    }

    fn group_params<T: Float>(params: [T; 49]) -> LayerParams<T, 20> {
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

// Extract and validate the "MODEL" line of the file
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
fn parse_param<T: Float>(param: &str) -> Result<T, IOError> {
    if param.starts_with('-') {
        return Err(IOError::new(
            ErrorKind::InvalidData,
            format!("A negative parameter was found: {param}"),
        ));
    }
    let param = param.parse::<T>().map_err(|_| {
        IOError::new(
            ErrorKind::InvalidData,
            format!("When parsing the model parameters, the value {param} could not be parsed as a float"),
        )
    })?;
    let mut param = -param.ln();
    if param.is_nan() {
        param = T::INFINITY;
    }
    Ok(param)
}

struct LayerIter<'a, T, const S: usize, const L: usize> {
    raw_layers: RawLayerIter<'a, T, S, L>,
    last_layer: LayerParams<T, S>,
}

impl<'a, T: Float, const S: usize, const L: usize> LayerIter<'a, T, S, L>
where
    SupportedConfig: ParserConfig<S, L>,
{
    fn new(lines: &'a mut LineIterator) -> Result<Self, IOError> {
        let mut raw_layers = RawLayerIter::new(lines);
        let last_layer = raw_layers
            .next()
            .ok_or(IOError::new(ErrorKind::InvalidData, "No model layers found!"))??;
        Ok(Self { raw_layers, last_layer })
    }
}

impl<T: Float, const S: usize, const L: usize> Iterator for LayerIter<'_, T, S, L>
where
    SupportedConfig: ParserConfig<S, L>,
{
    type Item = Result<LayerParams<T, S>, IOError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut next_layer = unwrap_or_return_some_err!(self.raw_layers.next()?);

        for param in PARAMS_TO_COPY {
            self.last_layer.transition[param] = next_layer.transition[param];
        }
        self.last_layer.emission_match = std::mem::replace(&mut next_layer.emission_match, EmissionParams([T::ZERO; S]));

        Some(Ok(std::mem::replace(&mut self.last_layer, next_layer)))
    }
}

struct RawLayerIter<'a, T, const S: usize, const L: usize> {
    lines:          &'a mut LineIterator,
    finished:       bool,
    expected_layer: Option<usize>,
    phantom:        PhantomData<T>,
}

impl<'a, T, const S: usize, const L: usize> RawLayerIter<'a, T, S, L> {
    fn new(lines: &'a mut LineIterator) -> Self {
        Self {
            lines,
            finished: false,
            expected_layer: Some(0),
            phantom: PhantomData,
        }
    }
}

impl<T: Float, const S: usize, const L: usize> Iterator for RawLayerIter<'_, T, S, L>
where
    SupportedConfig: ParserConfig<S, L>,
{
    type Item = Result<LayerParams<T, S>, IOError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        let line = unwrap_or_return_some_err!(self.lines.next().transpose());

        let Some(line) = line else {
            return Some(Err(IOError::new(ErrorKind::InvalidData, "File ended before ENDMODEL")));
        };

        let mut tokens = line.split_whitespace();
        let token = tokens.next().unwrap_or("");

        if token.eq_ignore_ascii_case("ENDMODEL") {
            self.finished = true;
            return None;
        }

        let Some(expected_layer) = self.expected_layer else {
            return Some(Err(IOError::new(
                ErrorKind::InvalidData,
                format!("Unexpected token {token} found between END layer and ENDMODEL"),
            )));
        };

        if token.eq_ignore_ascii_case("FREQAVE") {
            unwrap_or_return_some_err!(SupportedConfig::parse_layer_params::<T>(tokens, self.lines));
            return self.next();
        }

        let layer_number = unwrap_or_return_some_err!(parse_layer_name(token));

        if let Some(layer_number) = layer_number {
            if layer_number != expected_layer {
                return Some(Err(IOError::new(
                    ErrorKind::InvalidData,
                    format!("Found model node {layer_number}, expected model node {expected_layer}"),
                )));
            }

            self.expected_layer = Some(expected_layer + 1);
        } else {
            self.expected_layer = None;
        }

        let layer = match SupportedConfig::parse_layer_params(tokens, self.lines) {
            Ok(layer) => layer,
            Err(e) => return Some(Err(e)),
        };

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
                self.finished = true;
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
