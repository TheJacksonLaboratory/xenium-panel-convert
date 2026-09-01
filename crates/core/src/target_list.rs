use std::collections::{HashMap, HashSet};

use chemistry::{EnsemblId, GeneName, UnvalidatedEnsemblId};

use crate::{
    error::Hinted,
    target_list::{
        csv_util::{read_csv_trimmed, rename_fields},
        error::{TargetError, TargetErrorSet},
        target::{UnvalidatedTarget, ValidTarget},
    },
};

pub mod chemistry;
mod csv_util;
pub mod error;
pub mod target;
pub mod xenium_panel_designer;

pub fn parse_target_list(
    target_list: &str,
    field_aliases: &HashMap<&str, &str>,
    ensembl_id_to_gene: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)> + Copy,
) -> Result<Vec<ValidTarget>, Vec<TargetErrorSet>> {
    const N_GENES: usize = 500;

    let mut reader = read_csv_trimmed(target_list);
    // If we can't get headers, just return early
    let headers = reader.headers().map_err(|e| {
        vec![TargetErrorSet {
            line_number: None,
            submitted_target: None,
            errors: vec![Hinted::new(TargetError::from(e))],
        }]
    })?;

    let fieldnames = rename_fields(headers, field_aliases);

    let mut valid_targets = Vec::with_capacity(N_GENES);
    let mut seen_genes = HashSet::with_capacity(N_GENES);
    let mut errors = Vec::new();

    for record in reader.records() {
        let record = match record {
            Ok(record) => record,
            Err(err) => {
                errors.push(TargetErrorSet {
                    line_number: None,
                    submitted_target: None,
                    errors: vec![Hinted::new(TargetError::from(err))],
                });

                continue;
            }
        };

        let line_number = record.position().map(csv::Position::line);
        let submitted_target = UnvalidatedTarget::from_record(&record, &fieldnames);

        let row_errors = match ValidTarget::from_unvalidated(&submitted_target, ensembl_id_to_gene)
        {
            Ok(valid_target) => {
                // Cloning is cheap for the vast majority of genes
                if seen_genes.insert(valid_target.gene().clone()) {
                    valid_targets.push(valid_target);

                    continue;
                }

                vec![TargetError::DuplicateGene]
            }
            Err(row_errors) => row_errors,
        };

        errors.push(TargetErrorSet {
            line_number,
            submitted_target: Some(submitted_target),
            errors: row_errors.into_iter().map(Hinted::new).collect(),
        });
    }

    if !errors.is_empty() {
        return Err(errors);
    }

    Ok(valid_targets)
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use crate::{
        error::Hinted,
        target_list::{
            TargetError,
            chemistry::{tests::tp53_ensembl_id, xenium_v1_human_ensembl_id_to_gene},
            parse_target_list,
        },
    };

    #[test]
    fn valid_target_list() {
        let gene_list = "ensembl_id,gene_name,group,priority\nENSG00000141510,TP53,group0,\
                         must_have\nENSG00000116678,LEPR,group0,desired\nENSG00000120802,TMPO,\
                         group1,backup";

        let targets = parse_target_list(
            gene_list,
            &HashMap::new(),
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap();

        let gene_names: Vec<_> = targets
            .iter()
            .map(|target| {
                let (_, gene_name) = target.gene().as_strs();
                gene_name.unwrap()
            })
            .collect();
        assert_eq!(
            gene_names,
            ["TP53", "LEPR", "TMPO"],
            "targets were not returned in the order they were submitted"
        );
    }

    #[test]
    fn duplicate_targets() {
        let ensembl_id = tp53_ensembl_id();
        let ensembl_id_str = ensembl_id.as_str();

        // Two idential rows. We have to split this into 3 lines because cargo +nightly
        // fmt destroys it otherwise
        let header = "ensembl_id,gene_name,group,priority";
        let row = format!("{ensembl_id_str},TP53,group0,must_have");
        let gene_list = format!("{header}\n{row}\n{row}");

        let errors = parse_target_list(
            &gene_list,
            &HashMap::new(),
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        assert_eq!(errors.len(), 1, "did not find exactly 1 error");
        assert_eq!(errors[0].errors, [Hinted::new(TargetError::DuplicateGene)]);
    }

    #[test]
    fn error_reports_correct_file_line_number() {
        let gene_list = "ensembl_id,gene_name,group,priority\nid,gene,0,must_have";

        let errors = parse_target_list(
            gene_list,
            &HashMap::new(),
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        assert_eq!(errors.len(), 1);
        assert_eq!(errors[0].line_number, Some(2));
    }
}
