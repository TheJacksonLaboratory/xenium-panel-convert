use std::{fs, path::PathBuf};

use serde::Deserialize;

use crate::common::{parse_gene_list_from_csv, write_map_to_file};

const XENIUM_V1_GENE_LIST: &str = "xenium-gene-lists/human_and_mouse_2020-A-ref-yesprobe-genes_v1assay.csv";

const XENIUM_PRIME_GENE_LIST: &str = "xenium-gene-lists/human_and_mouse_genes_xenium_prime-yesprobe-genes.csv";

pub(crate) fn write_gene_maps() -> anyhow::Result<()> {
    let gene_lists: Vec<_> = [XENIUM_V1_GENE_LIST, XENIUM_PRIME_GENE_LIST]
        .map(|p| parse_gene_list_from_csv(&fs::read(p).unwrap()))
        .into_iter()
        .collect::<Result<_, _>>()?;

    let gene_maps: Vec<_> = gene_lists.iter().map(|list| construct_maps(list)).collect();

    let [
        GeneMaps {
            homo_sapiens: v1_human,
            mus_musculus: v1_mouse,
        },
        GeneMaps {
            homo_sapiens: prime_human,
            mus_musculus: prime_mouse,
        },
    ] = gene_maps.as_array().unwrap();

    for (filename, map_name, gene_map) in [
        ("xenium_v1_human.rs", "XENIUM_V1_HUMAN_GENES", v1_human),
        ("xenium_v1_mouse.rs", "XENIUM_V1_MOUSE_GENES", v1_mouse),
        ("xenium_prime_human.rs", "XENIUM_PRIME_HUMAN_GENES", prime_human),
        ("xenium_prime_mouse.rs", "XENIUM_PRIME_MOUSE_GENES", prime_mouse),
    ] {
        write_map_to_file(
            &PathBuf::from(format!("src/target_list/chemistry/{filename}")),
            map_name,
            gene_map,
        )?;
    }

    Ok(())
}

#[derive(Deserialize, Clone)]
struct Gene {
    #[serde(rename = "Species")]
    species: String,
    #[serde(rename = "Ensembl ID")]
    ensembl_id: String,
    #[serde(rename = "Gene symbol")]
    symbol: String,
}

fn construct_maps(gene_list: &[Gene]) -> GeneMaps<'_> {
    fn insert_gene<'a>(ensembl_id: &'a str, gene_symbol: &'a str, map: &mut phf_codegen::Map<'a, &'a str>) {
        map.entry(ensembl_id, format!(r#""{gene_symbol}""#));
    }

    let mut homo_sapiens = phf_codegen::Map::new();
    let mut mus_musculus = phf_codegen::Map::new();

    for Gene {
        species,
        ensembl_id,
        symbol,
    } in gene_list
    {
        match species.as_str() {
            "Homo sapiens" => insert_gene(ensembl_id, symbol, &mut homo_sapiens),
            "Mus musculus" => insert_gene(ensembl_id, symbol, &mut mus_musculus),
            s => panic!("species {s} not expected"),
        }
    }

    GeneMaps {
        homo_sapiens,
        mus_musculus,
    }
}

struct GeneMaps<'a> {
    homo_sapiens: phf_codegen::Map<'a, &'a str>,
    mus_musculus: phf_codegen::Map<'a, &'a str>,
}
