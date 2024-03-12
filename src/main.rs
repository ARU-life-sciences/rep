use rep::cli::pipeline;
use rep::error::Result;

// we want to build an application which will execute
// the repeatmodeler pipeline, and custom scripts
// this will make much use of `Command` to use shell
// commands to run the pipeline

fn main() -> Result<()> {
    let args = match pipeline() {
        Ok(_) => (),
        Err(err) => {
            eprintln!("{err}");
            std::process::exit(1);
        }
    };

    println!("{:?}", args);

    Ok(())
}
