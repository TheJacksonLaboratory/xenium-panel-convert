#![expect(clippy::doc_markdown)]
use std::{collections::HashMap, fs};

use anyhow::{Context, anyhow};
use camino::{Utf8Path, Utf8PathBuf};
use xenium_panel_convert_core::target_list::{
    TargetList,
    chemistry::{
        Chemistry, Species, xenium_prime_human_ensembl_id_to_gene, xenium_prime_mouse_ensembl_id_to_gene,
        xenium_v1_human_ensembl_id_to_gene, xenium_v1_mouse_ensembl_id_to_gene,
    },
    parse_target_list,
    target::ValidTargetCsv,
    xenium_panel_designer::XeniumPanelDesignerCsv,
};

use crate::write::{write_csv_with_header_to_file, write_json_to_file};

pub(super) fn convert_target_list(
    TargetListCliOptions {
        targets_path,
        field_alias_file,
        field_aliases: field_aliases_from_cli,
        species,
        chemistry,
    }: &TargetListCliOptions,
    output_dir: &Utf8Path,
) -> anyhow::Result<Option<TargetList>> {
    let target_list =
        fs::read_to_string(targets_path).with_context(|| format!("failed to read target-list from {targets_path}"))?;

    let field_aliases_from_file = read_field_aliases_from_file(field_alias_file.as_deref())?;

    let field_aliases = combine_field_aliases(&field_aliases_from_file, field_aliases_from_cli)?;

    let ensembl_id_to_gene = match (species, chemistry) {
        (Species::HomoSapiens, Chemistry::V1) => xenium_v1_human_ensembl_id_to_gene,
        (Species::HomoSapiens, Chemistry::Prime) => xenium_prime_human_ensembl_id_to_gene,
        (Species::MusMusculus, Chemistry::V1) => xenium_v1_mouse_ensembl_id_to_gene,
        (Species::MusMusculus, Chemistry::Prime) => xenium_prime_mouse_ensembl_id_to_gene,
    };

    let result = parse_target_list(&target_list, &field_aliases, ensembl_id_to_gene);

    let output_file_path = |filename| output_dir.join(filename);

    match result {
        Ok(target_list) => {
            let valid_target_csv = ValidTargetCsv::from_target_list(&target_list);
            write_csv_with_header_to_file(
                valid_target_csv.header(),
                valid_target_csv.rows(),
                &output_file_path("validated-targets.csv"),
            )?;

            let xpd_csv = XeniumPanelDesignerCsv::from_valid_targets(target_list.targets());
            write_csv_with_header_to_file(
                xpd_csv.header(),
                xpd_csv.rows(),
                &output_file_path("xenium-panel-designer-targets.csv"),
            )?;

            Ok(Some(target_list))
        }
        Err(e) => {
            write_json_to_file(&e, &output_file_path("target-list-errors.json"))?;

            Ok(None)
        }
    }
}

#[derive(Debug, Clone, clap::Args)]
pub(super) struct TargetListCliOptions {
    /// A path to a CSV-formatted file of panel targets. See https://github.com/TheJacksonLaboratory/xenium-panel-convert for more details.
    #[clap(long, short)]
    targets_path: Utf8PathBuf,
    /// A path to a TOML file mapping field-aliases to canonical fieldnames. See https://github.com/TheJacksonLaboratory/xenium-panel-convert for more details.
    #[clap(long, short = 'f')]
    field_alias_file: Option<Utf8PathBuf>,
    /// One or more field-aliases specified as '<ALIAS>=<FIELD>'. This is useful for processing CSV-files without
    /// editing their header.
    #[clap(long, short = 'a', value_parser = parse_field_aliases)]
    field_aliases: Vec<(String, String)>,
    /// The species of your target-list, used in conjunction with --chemistry to determine the set of genes allowed by
    /// the panel designer.
    #[clap(long, short)]
    pub(super) species: Species,
    /// The desired chemistry of your eventual panel, used in conjunction with --species to determine the set of genes
    /// allowed by the panel designer.
    #[clap(long, short)]
    chemistry: Chemistry,
}

fn parse_field_aliases(s: &str) -> anyhow::Result<(String, String)> {
    s.split_once('=')
        .map(|(alias, field)| (alias.to_owned(), field.to_owned()))
        .ok_or_else(|| anyhow!("field aliases must be specified as '<ALIAS>=<FIELD>'"))
}

fn read_field_aliases_from_file(path: Option<&Utf8Path>) -> anyhow::Result<Vec<u8>> {
    let Some(path) = path else {
        return Ok(Vec::new());
    };

    fs::read(path).with_context(|| format!("failed to read file {path}"))
}

fn combine_field_aliases<'a>(
    field_alias_file_contents: &'a [u8],
    field_aliases_from_cli: &'a [(String, String)],
) -> anyhow::Result<HashMap<&'a str, &'a str>> {
    let mut field_aliases: HashMap<&str, &str> = toml::from_slice(field_alias_file_contents)?;

    for (alias, field) in field_aliases_from_cli {
        field_aliases.insert(alias, field);
    }

    Ok(field_aliases)
}

#[cfg(test)]
mod tests {
    use crate::targets::{combine_field_aliases, parse_field_aliases};

    #[test]
    fn field_aliases_must_use_equals() {
        assert_eq!(
            parse_field_aliases("alias=field").unwrap(),
            ("alias".to_owned(), "field".to_owned())
        );

        assert!(parse_field_aliases("alias").is_err());
    }

    #[test]
    fn field_aliases_are_combined_correctly() {
        let field_aliases = ["alias1", "field1", "alias2", "field2"];

        let field_aliases: Vec<(String, String)> = field_aliases
            .chunks(2)
            .map(|alias_field| (alias_field[0].to_owned(), alias_field[1].to_owned()))
            .collect();

        let field_aliases = combine_field_aliases(br#"alias1 = "field2""#, &field_aliases).unwrap();

        assert_eq!(
            field_aliases,
            [("alias1", "field1"), ("alias2", "field2")].into_iter().collect()
        );
    }
}
