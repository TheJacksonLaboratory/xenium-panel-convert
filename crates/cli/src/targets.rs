#![allow(clippy::doc_markdown)]
use std::{collections::HashMap, fs};

use anyhow::{Context, anyhow};
use camino::{Utf8Path, Utf8PathBuf};
use xenium_panel_convert_core::target_list::{
    chemistry::{
        Chemistry, Species, xenium_prime_human_ensembl_id_to_gene,
        xenium_prime_mouse_ensembl_id_to_gene, xenium_v1_human_ensembl_id_to_gene,
        xenium_v1_mouse_ensembl_id_to_gene,
    },
    parse_target_list,
    target::{ValidTarget, to_valid_target_csv_rows},
    xenium_panel_designer::to_xenium_panel_designer_csv_rows,
};

use crate::write::{write_csv_to_file, write_json_to_file};

pub(super) fn convert_target_list(
    TargetListCliOptions {
        targets_path,
        field_alias_file,
        field_aliases: field_aliases_from_cli,
        species,
        chemistry,
    }: &TargetListCliOptions,
    output_dir: &Utf8Path,
) -> anyhow::Result<Option<Vec<ValidTarget>>> {
    let target_list = fs::read_to_string(targets_path)
        .with_context(|| format!("failed to read target-list from {targets_path}"))?;

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
        Ok(targets) => {
            write_csv_to_file(
                &to_valid_target_csv_rows(&targets),
                &output_file_path("validated-targets.csv"),
            )?;

            write_csv_to_file(
                &to_xenium_panel_designer_csv_rows(&targets),
                &output_file_path("xenium-panel-designer-targets.csv"),
            )?;

            Ok(Some(targets))
        }
        Err(e) => {
            write_json_to_file(&e, &output_file_path("target-list-errors.json"))?;

            Ok(None)
        }
    }
}

#[derive(Debug, Clone, clap::Args)]
pub(super) struct TargetListCliOptions {
    /// A path to a CSV-formatted file with the header
    /// "ensembl_id,gene_name,group,priority", where "priority" is one of
    /// "must_have", "desired", or "backup".
    #[clap(long, short)]
    targets_path: Utf8PathBuf,
    /// A path to a TOML file containing a map of field-aliases for the
    /// target-list CSV. The file should map aliases to canonical fieldnames.
    #[clap(long, short = 'f')]
    field_alias_file: Option<Utf8PathBuf>,
    /// One or more field-aliases specified as '<ALIAS>=<FIELD>'
    #[clap(long, short = 'a', value_parser = parse_field_aliases)]
    field_aliases: Vec<(String, String)>,
    #[clap(long, short)]
    pub(super) species: Species,
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
            [("alias1", "field1"), ("alias2", "field2")]
                .into_iter()
                .collect()
        );
    }
}
