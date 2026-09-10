use std::str::FromStr;

use crate::{complete_feature_sets::write_transcriptomes, xenium_panel_allowed_genes::write_gene_maps};

#[path = "build/common.rs"]
mod common;
#[path = "build/complete_feature_sets.rs"]
mod complete_feature_sets;
#[path = "build/xenium_panel_allowed_genes.rs"]
mod xenium_panel_allowed_genes;

fn main() -> anyhow::Result<()> {
    if is_true(option_env!("BUILD_XP_CONVERT_GENE_MAPS"))? {
        write_gene_maps()?;
    }

    if is_true(option_env!("BUILD_XP_CONVERT_TRANSCRIPTOMES"))? {
        write_transcriptomes()?;
    }

    Ok(())
}

fn is_true(env_var: Option<&'static str>) -> anyhow::Result<bool> {
    Ok(env_var.map(bool::from_str).transpose()?.is_some_and(|v| v))
}
