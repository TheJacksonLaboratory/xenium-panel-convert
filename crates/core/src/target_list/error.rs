use serde::Serialize;

use crate::{
    error::Hinted,
    target_list::{
        chemistry::{EnsemblId, GeneName},
        target::{UnvalidatedTarget, ValidGene},
    },
};

#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct TargetErrorSet {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub line_number: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub submitted_target: Option<UnvalidatedTarget>,
    pub errors: Vec<Hinted<TargetError>>,
}

#[derive(Clone, Debug, Serialize, PartialEq, thiserror::Error)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum TargetError {
    #[error("the CSV could not be parsed ({reason}) - ensure it is properly formatted")]
    MalformedCsv { reason: String },
    #[error("the field {fieldname} is missing - add it to the CSV")]
    MissingField { fieldname: &'static str },
    #[error("{value} is not a valid {field} - change it to one of {}", allowed.join(", "))]
    InvalidValue {
        field: &'static str,
        value: String,
        allowed: &'static [&'static str],
    },
    #[error(
        "the Ensembl ID is versioned or lowercase - remove the version and uppercase the ID"
    )]
    VersionedOrLowercaseEnsemblId { correct_gene: Option<ValidGene> },
    #[error("no Ensembl ID was provided - add one")]
    NoEnsemblId,
    #[error(
        "no gene name was provided - add one (based on the Ensembl ID, it is probably \
         {probable_gene_name})"
    )]
    NoGeneName { probable_gene_name: GeneName },
    #[error(
        "the gene name corresponding to the Ensembl ID {ensembl_id} is {correct_gene_name} - \
         change either the Ensembl ID or the gene name so they match"
    )]
    EnsemblIdGeneNameMismatch {
        ensembl_id: EnsemblId,
        correct_gene_name: GeneName,
    },
    #[error(
        "this gene is not available for the chosen chemistry - see the 10x Genomics allowed genes at: https://www.10xgenomics.com/support/software/xenium-panel-designer/latest/tutorials/create-gene-list#yesprobe"
    )]
    GeneNotFound,
    #[error("this gene appears more than once in the target-list - remove this entry")]
    DuplicateGene,
}

impl From<csv::Error> for TargetError {
    fn from(err: csv::Error) -> Self {
        Self::from(&err)
    }
}

impl<'a> From<&'a csv::Error> for TargetError {
    fn from(err: &'a csv::Error) -> Self {
        Self::MalformedCsv {
            reason: err.to_string(),
        }
    }
}
