use std::{
    error::Error as StdError,
    fmt,
    io::Error as IOError,
    io::ErrorKind as IOErrorKind,
    num::{ParseFloatError, ParseIntError},
};

use anyhow::Error as AnyhowError;
use csv::Error as CsvError;

// A type alias for `Result<T, rep::Error>`.
pub type Result<T> = anyhow::Result<T, Error>;

// An error that can happen.
#[derive(Debug)]
pub struct Error(Box<ErrorKind>);

impl Error {
    // A crate private constructor for `Error`.
    pub(crate) fn new(kind: ErrorKind) -> Error {
        Error(Box::new(kind))
    }

    // Return the specific type of this error.
    pub fn kind(&self) -> &ErrorKind {
        &self.0
    }

    // Unwrap this error into its underlying type.
    pub fn into_kind(self) -> ErrorKind {
        *self.0
    }
}

// The specific type of error that can occur.
#[derive(Debug)]
pub enum ErrorKind {
    IO(IOErrorKind),
    GenericCli(String),
    FastaReading(AnyhowError),
    Parsef64(ParseFloatError),
    ParseInt(ParseIntError),
    BlastParse(CsvError),
}

impl StdError for Error {}

// the display of the error. We can make this more fancy later
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &*self.0 {
            ErrorKind::IO(err) => write!(f, "IO error: {}", err),
            ErrorKind::GenericCli(msg) => write!(f, "Generic CLI error: {}", msg),
            ErrorKind::FastaReading(err) => write!(f, "Fasta reading error: {}", err),
            ErrorKind::Parsef64(err) => write!(f, "Error parsing float: {}", err),
            ErrorKind::ParseInt(err) => write!(f, "Error parsing int: {}", err),
            ErrorKind::BlastParse(err) => write!(f, "Error parsing BLAST output: {}", err),
        }
    }
}

// down here are all the conversion implementations
impl From<IOError> for Error {
    fn from(err: IOError) -> Error {
        Error::new(ErrorKind::IO(err.kind()))
    }
}

impl From<AnyhowError> for Error {
    fn from(err: AnyhowError) -> Error {
        Error::new(ErrorKind::FastaReading(err))
    }
}

impl From<ParseFloatError> for Error {
    fn from(err: ParseFloatError) -> Error {
        Error::new(ErrorKind::Parsef64(err))
    }
}

impl From<ParseIntError> for Error {
    fn from(err: ParseIntError) -> Error {
        Error::new(ErrorKind::ParseInt(err))
    }
}

impl From<CsvError> for Error {
    fn from(err: CsvError) -> Error {
        Error::new(ErrorKind::BlastParse(err))
    }
}
