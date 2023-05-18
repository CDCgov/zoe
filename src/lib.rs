#![warn(clippy::all, clippy::pedantic)]
#![allow(
    clippy::module_name_repetitions,
    clippy::similar_names,
    clippy::wildcard_imports,
    clippy::enum_glob_use
)]
#![feature(test, portable_simd, const_fn_floating_point_arithmetic)]
pub mod data;
pub mod distance;

pub(crate) mod math;
