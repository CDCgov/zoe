use crate::{
    alignment::phmm::{
        EmissionParams, LayerParams, Phmm,
        PhmmState::{self, *},
        TransitionParams,
    },
    data::mappings::UNAMBIG_DNA_PROFILE_MAP,
    math::Float,
};
use std::{
    fs::File,
    io::{BufRead, BufReader, Error as IOError, ErrorKind, Lines},
    path::Path,
};

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

impl<T: Float> Phmm<T, 4> {
    /// Parse a SAM model file representing a pHMM for a DNA alphabet.
    ///
    /// ## Errors
    /// The file must meet the specifications for a SAM model file.
    /// Additionally, it must be a model file (not a regularizer or null model)
    /// and it must have the alphabet specified as `DNA`.
    pub fn from_sam_file<P>(filename: P) -> Result<Self, std::io::Error>
    where
        P: AsRef<Path>, {
        let mut lines = BufReader::new(File::open(filename)?).lines();

        let mut line = get_next_line(&mut lines, "Could not locate initial line in file")?;
        // The unwrap will always be successful, since `get_next_line` returns a
        // non-empty trimmed string
        match line.split_whitespace().next().unwrap_or("") {
            "MODEL" => line = get_next_line(&mut lines, "Could not locate the model parameters in the file")?,
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

        let mapping = match line.strip_prefix("alphabet ") {
            Some("DNA") => &UNAMBIG_DNA_PROFILE_MAP,
            Some(_) => return Err(IOError::new(ErrorKind::InvalidData, "Unsupported alphabet specified")),
            None => return Err(IOError::new(ErrorKind::InvalidData, "Expected alphabet line")),
        };

        let mut layers = Vec::new();
        let mut current_layer = 0;
        let mut layers_left = None;
        let Some(mut last_layer) = get_next_dna_params::<T>(&mut lines, &mut current_layer, &mut layers_left)? else {
            return Err(IOError::new(
                ErrorKind::InvalidData,
                "At least one layer must be present in the model",
            ));
        };

        while let Some(mut layer) = get_next_dna_params::<T>(&mut lines, &mut current_layer, &mut layers_left)? {
            for param in PARAMS_TO_COPY {
                last_layer.transition[param] = layer.transition[param];
            }
            last_layer.emission_match = std::mem::replace(&mut layer.emission_match, EmissionParams([T::ZERO; 4]));
            layers.push(last_layer);
            last_layer = layer;
        }

        Ok(Phmm { mapping, params: layers })
    }
}

fn get_next_dna_params<T: Float>(
    lines: &mut Lines<BufReader<File>>, current_layer: &mut usize, layers_left: &mut Option<usize>,
) -> Result<Option<LayerParams<T, 4>>, std::io::Error> {
    let line = get_next_line(lines, "File ended before ENDMODEL")?;
    let mut tokens = line.split_whitespace();
    let Some(token) = tokens.next() else {
        return get_next_dna_params::<T>(lines, current_layer, layers_left);
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
        parse_dna_layer_params::<T>(tokens, lines)?;
        return get_next_dna_params(lines, current_layer, layers_left);
    }

    check_layer_name(token, current_layer, layers_left)?;
    parse_dna_layer_params(tokens, lines).map(Some)
}

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

// Check that the layer name is valid, given the expectation of the layer about
// to be parsed (`current_layer`) and the expectation for the number of layers
// left including the current one (`layers_left`). The function returns the
// updated values for `layers_left`.
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

/// Parses the parameters for a DNA pHMM layer from a set of tokens
/// (`rest_of_line`) followed by as many lines from `lines` as needed. There are
/// assumed to be 17 parameters, representing the 9 transition probabilities, 4
/// match emission probabilities, and 4 insert emission probabilities.
///
/// ## Errors
/// * `rest_of_line` and `lines` must contain sufficiently many tokens to fill a
///   full set of parameters
/// * All parameters must parse successfully, and there should be no extra
///   parameters on any line (see [`fill_params_from_iter`])
fn parse_dna_layer_params<'a, T: Float>(
    mut rest_of_line: impl Iterator<Item = &'a str>, lines: &mut Lines<BufReader<File>>,
) -> Result<LayerParams<T, 4>, std::io::Error> {
    let mut params = vec![T::ZERO; 17];

    let mut i = fill_params_from_iter::<T, 4>(rest_of_line.by_ref(), &mut params, 0)?;

    for line in lines {
        let line = line?;
        let mut split = line.split_whitespace();
        i = fill_params_from_iter::<T, 4>(split.by_ref(), &mut params, i)?;
        if i == 17 {
            break;
        }
    }

    if i < params.len() {
        return Err(IOError::new(
            ErrorKind::InvalidData,
            "The file ended before a full set of parameters was parsed",
        ));
    }

    let a = TransitionParams([
        [params[0], params[1], params[2]],
        [params[3], params[4], params[5]],
        [params[6], params[7], params[8]],
    ]);

    let e_m = EmissionParams([params[9], params[11], params[10], params[12]]);
    let e_i = EmissionParams([params[13], params[15], params[14], params[16]]);

    Ok(LayerParams {
        transition:      a,
        emission_match:  e_m,
        emission_insert: e_i,
    })
}

/// Uses the parameters in `iter` to fill spots in `params`, starting with `i`
/// and continuing until the end of `params` or `iter`. Returns the index of the
/// next spot to fill (or `params.len()` if it is full).
///
/// `params` is assumed to be of size $9+2S$.
///
/// ## Errors
/// * The parameters must successfully parse (see [`parse_param`])
/// * If `params` gets filled, then `iter` must not contain any more elements
fn fill_params_from_iter<'a, T: Float, const S: usize>(
    iter: &mut impl Iterator<Item = &'a str>, params: &mut [T], mut i: usize,
) -> Result<usize, std::io::Error> {
    for param in iter.take(9 + 2 * S - i) {
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

/// Parse a parameter from a string, and convert it to its negative natural
/// logarithm. A parameter that is 0 is converted to $\infty$.
///
/// ## Errors
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
