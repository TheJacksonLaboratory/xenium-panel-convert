#![expect(clippy::similar_names)]
use std::{
    collections::HashSet,
    fmt::Display,
    fs,
    path::{Path, PathBuf},
    str::FromStr,
};

use anyhow::Context;
use hdf5_metno::types::FixedAscii;
use serde::{Deserialize, Deserializer};

use crate::common::{parse_gene_list_from_csv, write_map_to_file};

pub(crate) fn write_transcriptomes() -> anyhow::Result<()> {
    write_3p_gene_lists()?;
    write_flex_gene_lists()
}

#[derive(Deserialize, Clone)]
struct Gene {
    #[serde(rename = "gene_id")]
    ensembl_id: String,
    probe_id: String,
    #[serde(deserialize_with = "deserialize_bool")]
    included: bool,
}

fn deserialize_bool<'de, D>(deserializer: D) -> Result<bool, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;

    Ok(bool::from_str(&s.to_lowercase()).expect("the string should always be TRUE or FALSE"))
}

fn write_3p_gene_lists() -> anyhow::Result<()> {
    let grch38_2020_a = (
        "datasets/SC3pv3_GEX_Human_PBMC_filtered_feature_bc_matrix.h5",
        36_601,
        "GRCH38_2020_A",
        "grch38_2020_a.rs",
    );

    let grch38_2024_a = (
        "datasets/5k_Human_Donor2_PBMC_3p_gem-x_5k_Human_Donor2_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.\
         h5",
        38_606,
        "GRCH38_2024_A",
        "grch38_2024_a.rs",
    );

    let mm10_2020_a = (
        "datasets/1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix.h5",
        32_285,
        "MM10_2020_A",
        "mm10_2020_a.rs",
    );

    let grcm39_2024_a = (
        "datasets/SOD1_G93A_mouse_spinal_cord_P112_specimen_1_SOD1_G93A_mouse_spinal_cord_P112_specimen_1_sample_filtered_feature_bc_matrix.h5",
        33_696,
        "GRCM39_2024_A",
        "grcm39_2024_a.rs",
    );

    for (h5_path, expected_n_genes, map_name, map_path) in [grch38_2020_a, grch38_2024_a, mm10_2020_a, grcm39_2024_a] {
        let (ensembl_ids, gene_symbols) = read_transcriptome_from_h5(h5_path)?;

        let map = construct_map(ensembl_ids.iter().zip(&gene_symbols), expected_n_genes);

        write_map_to_file(&create_map_path(map_path), map_name, &map)?;
    }

    Ok(())
}

fn write_flex_gene_lists() -> anyhow::Result<()> {
    // https://www.10xgenomics.com/support/flex-gene-expression/documentation/steps/probe-sets/chromium-frp-human-transcriptome-probe-set
    let grch38_2020_a_flex = (
        "probe-sets/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv",
        18_082,
        "GRCH38_2020_A_FLEX",
        "grch38_2020_a_flex.rs",
    );

    // https://www.10xgenomics.com/support/flex-gene-expression/documentation/steps/probe-sets/chromium-frp-human-transcriptome-probe-set-1-1
    let grch38_2024_a_flex_1_1 = (
        "probe-sets/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv",
        18_129,
        "GRCH38_2024_A_FLEX_V1_1",
        "grch38_2024_a_flex_v1_1.rs",
    );

    // https://www.10xgenomics.com/support/flex-gene-expression/documentation/steps/probe-sets/chromium-frp-human-transcriptome-probe-set-2-0
    let grch38_2024_a_flex_2_0 = (
        "probe-sets/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv",
        18_132,
        "GRCH38_2024_A_FLEX_V2_0",
        "grch38_2024_a_flex_v2_0.rs",
    );

    // https://www.10xgenomics.com/support/flex-gene-expression/documentation/steps/probe-sets/chromium-frp-mouse-transcriptome-probe-set
    let mm10_2020_a_flex = (
        "probe-sets/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv",
        19_059,
        "MM10_2020_A_FLEX",
        "mm10_2020_a_flex.rs",
    );

    // https://www.10xgenomics.com/support/flex-gene-expression/documentation/steps/probe-sets/chromium-frp-mouse-transcriptome-probe-set-1-1
    let grcm39_2024_a_flex_1_1 = (
        "probe-sets/Chromium_Mouse_Transcriptome_Probe_Set_v1.1.1_GRCm39-2024-A.csv",
        19_070,
        "GRCM39_2024_A_FLEX_V1_1",
        "grcm39_2024_a_flex_v1_1.rs",
    );

    // https://www.10xgenomics.com/support/flex-gene-expression/documentation/steps/probe-sets/chromium-frp-mouse-transcriptome-probe-set-2-0
    let grcm39_2024_a_flex_2_0 = (
        "probe-sets/Chromium_Mouse_Transcriptome_Probe_Set_v2.0.0_GRCm39-2024-A.csv",
        19_070,
        "GRCM39_2024_A_FLEX_V2_0",
        "grcm39_2024_a_flex_v2_0.rs",
    );

    for (raw_gene_list_path, expected_n_genes, map_name, map_path) in [
        grch38_2020_a_flex,
        grch38_2024_a_flex_1_1,
        grch38_2024_a_flex_2_0,
        mm10_2020_a_flex,
        grcm39_2024_a_flex_1_1,
        grcm39_2024_a_flex_2_0,
    ] {
        let genes: Vec<Gene> = parse_gene_list_from_csv(&fs::read(raw_gene_list_path)?)
            .with_context(|| format!("failed to parse gene-list for {map_path}"))?;
        let genes: Vec<_> = genes
            .iter()
            .filter(|g| g.included)
            .map(|g| {
                (
                    // These clones are necessary for reasons
                    g.ensembl_id.clone(),
                    g.probe_id.split('|').nth(1).map(str::to_owned).unwrap(),
                )
            })
            .collect();

        let gene_map = construct_map(genes.iter().map(|(e, g)| (e, g)), expected_n_genes);

        let map_path = format!("src/reference_dataset/transcriptome/{map_path}");
        write_map_to_file(Path::new(&map_path), map_name, &gene_map)?;
    }

    Ok(())
}

fn create_map_path(filename: &str) -> PathBuf {
    PathBuf::from(format!("src/reference_dataset/transcriptome/{filename}"))
}

// Human Ensembl IDs are 15 characters while mouse Ensembl IDs are 18
type EnsemblId = FixedAscii<18>;

// No gene name is likely to exceed 32 characters
type GeneName = FixedAscii<32>;

fn read_transcriptome_from_h5(file_path: &str) -> anyhow::Result<(Vec<EnsemblId>, Vec<GeneName>)> {
    let file = hdf5_metno::File::open(file_path)?;

    let ensembl_ids = file.dataset("matrix/features/id").and_then(|ds| ds.read_raw())?;
    let gene_symbols = file.dataset("matrix/features/name").and_then(|ds| ds.read_raw())?;

    Ok((ensembl_ids, gene_symbols))
}

fn construct_map<'a, EnsemblId, GeneName>(
    genes: impl Iterator<Item = (&'a EnsemblId, &'a GeneName)>,
    expected_n_genes: usize,
) -> phf_codegen::Map<'a, &'a str>
where
    EnsemblId: AsRef<str> + 'a,
    GeneName: Display + 'a,
{
    let mut seen = HashSet::with_capacity(65_536);
    let mut map = phf_codegen::Map::new();

    for (id, name) in genes {
        if seen.insert(id.as_ref()) {
            map.entry(id.as_ref(), format!(r#""{name}""#));
        }
    }

    assert_eq!(seen.len(), expected_n_genes, "wrong number of genes");

    map
}
