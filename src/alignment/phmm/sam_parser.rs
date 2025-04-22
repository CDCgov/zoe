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
    /// The file must meet the specifications for a SAM model file.
    /// Additionally, it must be a model file (not a regularizer or null model)
    /// and the specified alphabet must be dna.
    pub fn parse_dna_model(filename: impl AsRef<Path>) -> Result<Phmm<f32, 4>, std::io::Error> {
        SupportedConfig::parse_sam_file(filename)
    }

    /// Parse a protein pHMM from a SAM model file.
    ///
    /// ## Errors
    ///
    /// The file must meet the specifications for a SAM model file.
    /// Additionally, it must be a model file (not a regularizer or null model)
    /// and the specified alphabet must be protein.
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
    /// The file must meet the specifications for a SAM model file.
    /// Additionally, it must be a model file (not a regularizer or null model)
    /// and the specified alphabet must be dna.
    fn parse_dna_model(filename: impl AsRef<Path>) -> Result<Phmm<T, 4>, std::io::Error> {
        SupportedConfig::parse_sam_file(filename)
    }

    /// Parse a protein pHMM from a SAM model file.
    ///
    /// ## Errors
    ///
    /// The file must meet the specifications for a SAM model file.
    /// Additionally, it must be a model file (not a regularizer or null model)
    /// and the specified alphabet must be protein.
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
    fn get_mapping(mapping: &str) -> Result<&'static ByteIndexMap<S>, std::io::Error>;

    /// Given a flat array of `L=2*S+9` parameters, converts them into
    /// [`LayerParams`].
    fn group_params<T: Float>(params: [T; L]) -> LayerParams<T, S>;

    /// Parse a SAM model file representing a pHMM.
    ///
    /// ## Errors
    ///
    /// The file must meet the specifications for a SAM model file.
    /// Additionally, it must be a model file (not a regularizer or null model)
    /// and the specified alphabet must be compatible with the selected `S` and
    /// `L`.
    fn parse_sam_file<T: Float, P>(filename: P) -> Result<Phmm<T, S>, std::io::Error>
    where
        P: AsRef<Path>,
        SupportedConfig: ParserConfig<S, L>, {
        let mut lines = BufReader::new(File::open(filename)?).lines();

        let mut line = get_next_line(&mut lines, "Could not locate initial line in file")?;
        // The unwrap will always be successful, since `get_next_line` returns a
        // non-empty trimmed string
        match line.split_whitespace().next().unwrap_or("") {
            "MODEL" => {}
            "REGULARIZER" => {
                return Err(IOError::new(
                    ErrorKind::InvalidData,
                    "REGULARIZER is not supported by this parser",
                ));
            }
            "NULLMODEL" => {
                return Err(IOError::new(
                    ErrorKind::InvalidData,
                    "NULLMODEL is not supported by this parser",
                ));
            }
            _ => {
                return Err(IOError::new(ErrorKind::InvalidData, "Could not locate initial line in file"));
            }
        }

        line = get_next_line(&mut lines, "Could not locate alphabet line in file")?;
        line.make_ascii_uppercase();

        let Some(mut alphabet) = line.strip_prefix("ALPHABET") else {
            return Err(IOError::new(ErrorKind::InvalidData, "Expected alphabet line"));
        };
        alphabet = alphabet.trim();
        let mapping = SupportedConfig::get_mapping(alphabet)?;

        let mut layers = Vec::new();
        let mut current_layer = 0;
        let mut layers_left = None;
        let Some(mut last_layer) = SupportedConfig::get_next_layer::<T>(&mut lines, &mut current_layer, &mut layers_left)?
        else {
            return Err(IOError::new(
                ErrorKind::InvalidData,
                "At least one layer must be present in the model",
            ));
        };

        while let Some(mut layer) = SupportedConfig::get_next_layer::<T>(&mut lines, &mut current_layer, &mut layers_left)? {
            for param in PARAMS_TO_COPY {
                last_layer.transition[param] = layer.transition[param];
            }
            last_layer.emission_match = std::mem::replace(&mut layer.emission_match, EmissionParams([T::ZERO; S]));
            layers.push(last_layer);
            last_layer = layer;
        }

        Ok(Phmm { mapping, params: layers })
    }

    /// Retrieve the parameters for the next layer, ensuring that the layer is
    /// properly named and skipping `FREQAVE`. Once `ENDMODEL` is reached, this
    /// function returns `None`.
    fn get_next_layer<T: Float>(
        lines: &mut Lines<BufReader<File>>, current_layer: &mut usize, layers_left: &mut Option<usize>,
    ) -> Result<Option<LayerParams<T, S>>, std::io::Error> {
        let line = get_next_line(lines, "File ended before ENDMODEL")?;
        let mut tokens = line.split_whitespace();
        let Some(token) = tokens.next() else {
            return Self::get_next_layer::<T>(lines, current_layer, layers_left);
        };
        if token == "ENDMODEL" {
            if let Some(layers_left) = layers_left
                && *layers_left > 0
            {
                return Err(IOError::new(
                    ErrorKind::InvalidData,
                    "Some layers were missing from the end of the file",
                ));
            }
            return Ok(None);
        }

        if token == "FREQAVE" {
            Self::parse_layer_params::<T>(tokens, lines)?;
            return Self::get_next_layer(lines, current_layer, layers_left);
        }

        check_layer_name(token, current_layer, layers_left)?;
        Self::parse_layer_params(tokens, lines).map(Some)
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
        mut rest_of_line: impl Iterator<Item = &'a str>, lines: &mut Lines<BufReader<File>>,
    ) -> Result<LayerParams<T, S>, std::io::Error> {
        let mut params = [T::ZERO; L];

        let mut i = Self::fill_params_from_iter::<T>(rest_of_line.by_ref(), &mut params, 0)?;

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
    fn get_mapping(mapping: &str) -> Result<&'static ByteIndexMap<4>, std::io::Error> {
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
    fn get_mapping(mapping: &str) -> Result<&'static ByteIndexMap<20>, std::io::Error> {
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

/// Retrieve the next line from `lines`, giving a specified error message if
/// there are no more lines present. Leading and trailing whitespace is removed
/// from lines, and any empty lines or comment lines are skipped.
fn get_next_line(lines: &mut Lines<BufReader<File>>, message: &str) -> Result<String, std::io::Error> {
    loop {
        let line = lines
            .next()
            .transpose()?
            .ok_or(IOError::new(ErrorKind::InvalidData, message))?
            .trim()
            .to_string();

        if !line.is_empty() && !line.starts_with('%') {
            return Ok(line);
        }
    }
}

/// Check that the layer name is valid, given the expectation of the layer about
/// to be parsed (`current_layer`) and the expectation for the number of layers
/// left including the current one (`layers_left`). The function returns the
/// updated values for `layers_left`.
fn check_layer_name(token: &str, current_layer: &mut usize, layers_left: &mut Option<usize>) -> Result<(), IOError> {
    if token == "Begin" {
        if *current_layer == 0 {
            *current_layer += 1;
            return Ok(());
        }
    } else if token == "End" {
        if layers_left.is_none() || *layers_left == Some(1) {
            *current_layer += 1;
            *layers_left = Some(0);
            return Ok(());
        }
    } else if let Some(pos_token) = token.strip_prefix('-')
        && let Ok(layer) = pos_token.parse::<usize>()
    {
        if let Some(layers_left) = layers_left {
            if layer == *layers_left && *layers_left > 0 {
                *current_layer += 1;
                *layers_left -= 1;
                return Ok(());
            }
        } else {
            *current_layer += 1;
            *layers_left = Some(layer - 1);
            return Ok(());
        }
    } else if let Ok(layer) = token.parse::<usize>() {
        if let Some(layers_left) = layers_left {
            if *layers_left > 0 && *current_layer == layer {
                *current_layer += 1;
                *layers_left -= 1;
                return Ok(());
            }
        } else if *current_layer == layer {
            *current_layer += 1;
            return Ok(());
        }
    } else {
        return Err(IOError::new(ErrorKind::InvalidData, "Could not parse the layer name"));
    }

    Err(IOError::new(ErrorKind::InvalidData, "The layers were specified out of order"))
}

/// Parse a parameter from a string, and convert it to its negative natural
/// logarithm. A parameter that is 0 is converted to $\infty$.
///
/// ## Errors
///
/// * Negative parameters are not allowed
/// * The parameter must succeed when parsing to type `T`
fn parse_param<T: Float>(param: &str) -> Result<T, std::io::Error> {
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
