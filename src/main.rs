use rep::{pipeline, Result};

fn main() -> Result<()> {
    let result = match std::env::args().nth(1).as_deref() {
        Some("recover-copies") => {
            rep::cli::parse_recover_copies_args().and_then(|c| rep::run_copy_recovery(&c).map(drop))
        }
        Some("classify") => rep::cli::parse_classify_args()
            .and_then(|c| rep::classify::run_classify(&c))
            .map(|r| {
                eprintln!(
                    "Classified {} families: {:?}",
                    r.summary.classified_families, r.summary.by_disposition
                )
            }),
        Some("refdb") => rep::cli::run_refdb(),
        Some("discover") => rep::cli::run_discover(),
        Some("relatives") => rep::cli::run_relatives(),
        Some("relabel") => rep::cli::run_relabel(),
        Some("annotate") => rep::cli::run_annotate(),
        Some("curate") => rep::cli::run_curate(),
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
