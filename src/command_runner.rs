use crate::Result;
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
