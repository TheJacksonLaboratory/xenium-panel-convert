#![expect(clippy::doc_markdown)]
use std::{fs, process};

use anyhow::Context;
use camino::Utf8PathBuf;
use clap::Parser;

use crate::{
    all::convert_target_list_and_reference_datasets,
    reference_datasets::{ReferenceDatasetCliOptions, convert_all_reference_datasets_and_write, dataset_name},
    targets::{TargetListCliOptions, convert_target_list},
};

mod all;
mod reference_datasets;
mod targets;
mod write;

fn main() -> anyhow::Result<()> {
    let Cli { command, output_dir } = Cli::parse();

    fs::create_dir_all(&output_dir).with_context(|| format!("failed to create output directory {output_dir}"))?;

    match command {
        Command::Targets(options) => {
            convert_target_list(&options, &output_dir)?;
        }
        Command::References(options) => {
            convert_all_reference_datasets_and_write(&options, &output_dir)?;
        }
        Command::All {
            targets_options,
            references_options,
        } => convert_target_list_and_reference_datasets(&targets_options, &references_options, &output_dir)?,
        Command::Upgrade => {
            let updater_name = "xenium-panel-convert-update";

            let err_message1 = || format!("failed to run {updater_name}");

            let err_message2 = || {
                format!(
                    "failed to upgrade {} - go to {}/releases/latest to reinstall",
                    env!("CARGO_BIN_NAME"),
                    env!("CARGO_PKG_REPOSITORY")
                )
            };

            let output = process::Command::new(&updater_name)
                .output()
                .with_context(err_message1)
                .with_context(err_message2)?;

            let process::Output {
                status,
                mut stdout,
                stderr,
            } = output;

            if !stdout.is_empty() {
                stdout.push('\n' as u8);
            }

            let message = String::from_utf8(stdout.into_iter().chain(stderr).collect()).unwrap_or_default();

            if !status.success() {
                eprintln!("{message}");
            } else {
                println!("{message}");
            }
        }
    }

    Ok(())
}

#[derive(clap::Parser)]
#[clap(version)]
/// Convert files to the formats accepted by the 10x Genomics Xenium Panel Designer
struct Cli {
    #[clap(subcommand)]
    command: Command,
    /// The directory to write outputs into. If it doesn't exist, it is created.
    #[clap(long, short, global = true, default_value_t = Utf8PathBuf::from("xp-convert"))]
    output_dir: Utf8PathBuf,
}

#[derive(clap::Subcommand)]
enum Command {
    /// Convert a target-list to a format suitable for the Xenium Panel Designer.
    ///
    /// See https://github.com/TheJacksonLaboratory/xenium-panel-convert#xp-convert-targets for information on inputs and outputs.
    Targets(TargetListCliOptions),
    /// Convert scanpy-annotated single-cell RNA sequencing datasets to a format suitable for the Xenium Panel Designer.
    ///
    /// See https://github.com/TheJacksonLaboratory/xenium-panel-convert#xp-convert-references for information on inputs and outputs.
    References(ReferenceDatasetCliOptions),
    /// Convert a target-list and reference datasets to formats suitable for the Xenium Panel Designer.
    ///
    /// See https://github.com/TheJacksonLaboratory/xenium-panel-convert#xp-convert-all for information on inputs and outputs.
    All {
        #[clap(flatten)]
        targets_options: TargetListCliOptions,
        #[clap(flatten)]
        references_options: ReferenceDatasetCliOptions,
    },
    /// Update xp-convert to the latest version (runs xenium-panel-convert-update)
    Upgrade,
}
