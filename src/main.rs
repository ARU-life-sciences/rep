use rep::{pipeline, Result};

fn main() -> Result<()> {
    let result = match std::env::args().nth(1).as_deref() {
        Some("recover-copies") => {
            rep::cli::parse_recover_copies_args().and_then(|c| rep::run_copy_recovery(&c).map(drop))
        }
        _ => pipeline(),
    };
    match result {
        Ok(_) => (),
        Err(err) => {
            eprintln!("{err}");
            std::process::exit(1);
        }
    };
    Ok(())
}
