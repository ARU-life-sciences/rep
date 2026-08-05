use rep::{pipeline, Result};

fn main() -> Result<()> {
    if std::env::args().nth(1).as_deref() == Some("recover-copies") {
        let config = rep::cli::parse_recover_copies_args()?;
        rep::run_copy_recovery(&config)?;
        return Ok(());
    }
    match pipeline() {
        Ok(_) => (),
        Err(err) => {
            eprintln!("{err}");
            std::process::exit(1);
        }
    };
    Ok(())
}
