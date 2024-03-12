use rep::{pipeline, Result};

fn main() -> Result<()> {
    let args = match pipeline() {
        Ok(_) => (),
        Err(err) => {
            eprintln!("{err}");
            std::process::exit(1);
        }
    };

    Ok(())
}
