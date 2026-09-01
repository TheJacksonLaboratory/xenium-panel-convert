use serde::Serialize;

use crate::{
    common::ErrorVecExt,
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
    pub errors: Vec<TargetError>,
}

#[derive(Clone, Debug, Serialize, PartialEq, thiserror::Error)]
#[error("{error}")]
pub struct TargetError {
    #[serde(flatten)]
    pub error: TargetErrorInner,
    pub hint: String,
}

impl From<TargetErrorInner> for TargetError {
    fn from(error: TargetErrorInner) -> Self {
        Self {
            hint: error.to_string(),
            error,
        }
    }
}

#[derive(Clone, Debug, Serialize, PartialEq, thiserror::Error)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum TargetErrorInner {
    #[error("ensure the CSV is properly formatted")]
    MalformedCsv { reason: String },
    #[error("add the field {fieldname} to the CSV")]
    MissingField { fieldname: &'static str },
    #[error("change {value} to one of {}", allowed.join(","))]
    InvalidValue {
        field: &'static str,
        value: String,
        allowed: &'static [&'static str],
    },
    #[error("remove the Ensembl ID version and uppercase it")]
    VersionedOrLowercaseEnsemblId { correct_gene: Option<ValidGene> },
    #[error("add an Ensembl ID")]
    NoEnsemblId,
    #[error("add a gene name (based on the Ensembl ID, it is probably {probable_gene_name})")]
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
        "gene not found - see 10x Genomics allowed genes at: https://www.10xgenomics.com/support/software/xenium-panel-designer/latest/tutorials/create-gene-list#yesprobe"
    )]
    GeneNotFound,
    #[error("remove this entry from the gene-list")]
    DuplicateGene,
}

impl From<csv::Error> for TargetErrorInner {
    fn from(err: csv::Error) -> Self {
        Self::from(&err)
    }
}

impl<'a> From<&'a csv::Error> for TargetErrorInner {
    fn from(err: &'a csv::Error) -> Self {
        Self::MalformedCsv {
            reason: err.to_string(),
        }
    }
}

impl ErrorVecExt<TargetErrorInner> for Vec<TargetErrorInner> {
    fn push_err<T>(&mut self, err: TargetErrorInner) -> Option<T> {
        self.push(err);

        None
    }
}
