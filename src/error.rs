use std::{
    error::Error as StdError,
    fmt,
    io::Error as IOError,
    num::{ParseFloatError, ParseIntError},
};

use anyhow::Error as AnyhowError;
use csv::Error as CsvError;

// A type alias for `Result<T, rep::Error>`.
pub type Result<T> = std::result::Result<T, Error>;

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
    IO(IOError),
    GenericCli(String),
    Other(AnyhowError),
    Parsef64(ParseFloatError),
    ParseInt(ParseIntError),
    BlastParse(CsvError),
}

impl StdError for Error {
    fn source(&self) -> Option<&(dyn StdError + 'static)> {
        match &*self.0 {
            ErrorKind::IO(err) => Some(err),
            ErrorKind::GenericCli(_) => None,
            ErrorKind::Other(err) => Some(err.as_ref()),
            ErrorKind::Parsef64(err) => Some(err),
            ErrorKind::ParseInt(err) => Some(err),
            ErrorKind::BlastParse(err) => Some(err),
        }
    }
}

// the display of the error. We can make this more fancy later
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &*self.0 {
            ErrorKind::IO(err) => write!(f, "IO error: {}", err),
            ErrorKind::GenericCli(msg) => write!(f, "Generic CLI error: {}", msg),
            ErrorKind::Other(err) => write!(f, "Error: {}", err),
            ErrorKind::Parsef64(err) => write!(f, "Error parsing float: {}", err),
            ErrorKind::ParseInt(err) => write!(f, "Error parsing int: {}", err),
            ErrorKind::BlastParse(err) => write!(f, "Error parsing BLAST output: {}", err),
        }
    }
}

// down here are all the conversion implementations
impl From<IOError> for Error {
    fn from(err: IOError) -> Error {
        Error::new(ErrorKind::IO(err))
    }
}

impl From<AnyhowError> for Error {
    fn from(err: AnyhowError) -> Error {
        Error::new(ErrorKind::Other(err))
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
