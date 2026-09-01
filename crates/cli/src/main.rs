#![expect(clippy::doc_markdown)]
use std::fs;

use anyhow::Context;
use camino::Utf8PathBuf;
use clap::Parser;

use crate::{
    all::convert_target_list_and_reference_datasets,
    reference_datasets::{ReferenceDatasetCliOptions, convert_reference_datasets, dataset_name},
    targets::{TargetListCliOptions, convert_target_list},
};

mod all;
mod reference_datasets;
mod targets;
mod write;

fn main() -> anyhow::Result<()> {
    let Cli {
        command,
        output_dir,
    } = Cli::parse();

    fs::create_dir_all(&output_dir)
        .with_context(|| format!("failed to create output directory {output_dir}"))?;

    match command {
        Command::Targets(options) => {
            convert_target_list(&options, &output_dir)?;
        }
        Command::References(options) => {
            convert_reference_datasets(&options, &output_dir)?;
        }
        Command::All {
            targets_options,
            references_options,
        } => convert_target_list_and_reference_datasets(
            &targets_options,
            &references_options,
            &output_dir,
        )?,
    }

    Ok(())
}

#[derive(clap::Parser)]
#[clap(version)]
struct Cli {
    #[clap(subcommand)]
    command: Command,
    #[clap(long, short, global = true, default_value_t = Utf8PathBuf::from("xp-convert"))]
    output_dir: Utf8PathBuf,
}

#[derive(clap::Subcommand)]
enum Command {
    /// Convert a target-list to a format suitable for the Xenium Panel
    /// Designer.
    ///
    /// The target-list must be a CSV-file with the header:
    /// A "cleaned" version of the file will be saved at
    /// <OUTPUT_DIR>/validated-targets.csv, and the version for the panel
    /// designer will be saved at
    /// <OUTPUT_DIR>/xenium-panel-designer-target-list.csv. If errors are
    /// encountered, they are collected and written to
    /// <OUTPUT_DIR>/target-list-errors.json.
    Targets(TargetListCliOptions),
    /// Convert scanpy-annotated single-cell RNA sequencing datasets to a format
    /// suitable for the Xenium Panel Designer.
    ///
    /// Each annotated dataset (h5ad files generated with scanpy) is converted
    /// to a directory containing a matrix.h5 and annotations.csv. This
    /// directory can be fed to tar to create a compressed archive and uploaded
    /// directly to the panel designer as a reference dataset. If errors occur,
    /// they are collected and written to
    /// <OUTPUT_DIR>/<DATASET_PATH>-errors.json
    References(ReferenceDatasetCliOptions),
    /// Convert a target-list and reference datasets to formats suitable for the
    /// Xenium Panel Designer.
    ///
    /// This is the equivalent of running both commands at once with the added
    /// benefit of ensuring the genes in the target-list match the
    /// transcriptomes against which the reference datasets were aligned.
    All {
        #[clap(flatten)]
        targets_options: TargetListCliOptions,
        #[clap(flatten)]
        references_options: ReferenceDatasetCliOptions,
    },
}
