use std::{
    error::Error as StdError, fmt, io::Error as IOError, io::ErrorKind as IOErrorKind, result,
};

// A type alias for `Result<T, rep::Error>`.
pub type Result<T> = result::Result<T, Error>;

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
}

impl StdError for Error {}

// the display of the error. We can make this more fancy later
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &*self.0 {
            ErrorKind::IO(err) => write!(f, "IO error: {}", err),
            ErrorKind::GenericCli(msg) => write!(f, "Generic CLI error: {}", msg),
        }
    }
}

// down here are all the conversion implementations
impl From<IOError> for Error {
    fn from(err: IOError) -> Error {
        Error::new(ErrorKind::IO(err.kind()))
    }
}
