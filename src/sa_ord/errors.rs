use std::num::TryFromIntError;
use thiserror::Error;

/// Generic result type for the errors in this crate.
pub type Result<T, E = Error> = std::result::Result<T, E>;

/// Error types for this crate
#[derive(Error, Debug)]
pub enum Error {
    // General errors
    #[error("Suffix Array propagation tried to fetch a val ({val}) not in bin edges")]
    ValueNotInBinEdges { val: String },

    #[error(
        "Value marked as the beggining of an LMS string {beginning} didn't have corresponding end in mapping"
    )]
    LmsIndexError { beginning: usize },

    #[error(
        "The LMS substring starting at text index {text_index} was not found in mapping from substrings to reduced index"
    )]
    LmsReductionError { text_index: usize },

    #[error("Conversion error - {e}")]
    TryFromIntError { e: TryFromIntError },
}
