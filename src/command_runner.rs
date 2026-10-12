use crate::{Error, ErrorKind, Result};
use std::process::{Command, Output};

pub trait CommandRunner {
    fn run(&self, cmd: &mut Command) -> Result<Output>;
}

pub struct RealCommandRunner;

impl CommandRunner for RealCommandRunner {
    fn run(&self, cmd: &mut Command) -> Result<Output> {
        let output = cmd.output()?;
        Ok(output)
    }
}

pub(crate) fn run_checked(
    runner: &dyn CommandRunner,
    cmd: &mut Command,
    description: &str,
) -> Result<Output> {
    let output = runner.run(cmd).map_err(|e| match e.into_kind() {
        ErrorKind::IO(io) => Error::new(ErrorKind::GenericCli(format!(
            "{description}: could not run {}: {io}",
            cmd.get_program().to_string_lossy()
        ))),
        other => Error::new(other),
    })?;
    if !output.status.success() {
        return Err(Error::new(ErrorKind::GenericCli(format!(
            "{description} failed with status {}: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr).trim()
        ))));
    }
    Ok(output)
}
