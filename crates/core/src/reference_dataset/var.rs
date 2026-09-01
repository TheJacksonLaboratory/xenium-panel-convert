#![allow(clippy::result_large_err)]
use std::collections::HashSet;

use hdf5_metno::{
    File,
    types::{FixedAscii, VarLenUnicode},
};
use ndarray::Array1;
use serde::Serialize;

use crate::{
    error::collect_error,
    reference_dataset::{
        columns::{EnsemblIdCol, GeneNameCol},
        h5_util::{ReadH5FieldError, read_1d_string_dataset, to_ascii},
        transcriptome::Transcriptome,
    },
};

pub(super) fn read_features_from_h5ad(
    file: &File,
    ensembl_id_col: &EnsemblIdCol,
    gene_name_col: &GeneNameCol,
    transcriptome: Option<Transcriptome>,
) -> Result<Features, VarError> {
    let mut errors = Vec::new();

    let [Some(ensembl_ids), Some(gene_names), Some(feature_types)] = [
        ensembl_id_col.as_str(),
        gene_name_col.as_str(),
        "feature_types",
    ]
    .map(|column| read_1d_string_dataset(file, &format!("var/{column}")))
    .map(|field| collect_error(field, &mut errors)) else {
        return Err(VarError::InvalidH5Fields { errors });
    };

    check_feature_array_lens(&ensembl_ids, &gene_names, &feature_types)?;

    let features = Features {
        ensembl_ids: ensembl_ids.mapv(|s| to_ascii(&s)),
        gene_names: gene_names.mapv(|s| to_ascii(&s)),
        feature_types: feature_types.mapv(|s| to_ascii(&s)),
    };

    let Some(transcriptome) = transcriptome else {
        return Ok(features);
    };

    let n_genes_in_dataset = features.ensembl_ids.len();
    let expected_genes = transcriptome.gene_map(n_genes_in_dataset).ok_or_else(|| {
        let (n_expected_genes, n_expected_genes2) = transcriptome.n_genes();

        VarError::FilteredGenes {
            n_expected_genes,
            n_expected_genes2,
            n_found_genes: n_genes_in_dataset,
        }
    })?;

    validate_var_matches_transcriptome(&ensembl_ids, &gene_names, expected_genes)?;

    Ok(features)
}

fn check_feature_array_lens(
    ensembl_ids: &Array1<VarLenUnicode>,
    gene_names: &Array1<VarLenUnicode>,
    feature_types: &Array1<VarLenUnicode>,
) -> Result<(), VarError> {
    if ensembl_ids.len() != gene_names.len() || ensembl_ids.len() != feature_types.len() {
        return Err(VarError::InvalidShapes {
            ensembl_ids_len: ensembl_ids.len(),
            gene_names_len: gene_names.len(),
            feature_types_len: feature_types.len(),
        });
    }

    Ok(())
}

fn validate_var_matches_transcriptome(
    ensembl_ids: &Array1<VarLenUnicode>,
    gene_names: &Array1<VarLenUnicode>,
    expected_genes: &phf::Map<&str, &str>,
) -> Result<(), VarError> {
    let mut errors = Vec::new();
    let mut seen = HashSet::with_capacity(ensembl_ids.len());

    for (id, name) in ensembl_ids.iter().zip(gene_names) {
        if !seen.insert(id) {
            errors.push(VarRowError::DuplicateGene {
                ensembl_id: id.to_string(),
                gene_name: name.to_string(),
            });

            continue;
        }

        let Some(expected_gene_name) = expected_genes.get(id) else {
            errors.push(VarRowError::UnrecognizedEnsemblId {
                ensembl_id: id.to_string(),
            });

            continue;
        };

        if name != *expected_gene_name {
            errors.push(VarRowError::EnsemblIdGeneNameMismatch {
                ensembl_id: id.to_string(),
                expected_gene_name,
                found_gene_name: name.to_string(),
            });
        }
    }

    if errors.is_empty() {
        Ok(())
    } else {
        Err(VarError::Genes { errors })
    }
}

fn hints(errors: &[VarRowError]) -> Vec<&'static str> {
    let mut hints = Vec::new();

    for hint in errors.iter().map(VarRowError::hint) {
        if !hints.contains(&hint) {
            hints.push(hint);
        }
    }

    hints
}

// Human Ensembl IDs are 15 characters while mouse Ensembl IDs are 18
pub(super) type EnsemblId = FixedAscii<18>;

type EnsemblIds = Array1<EnsemblId>;

// No gene name is likely to exceed 32 characters
pub(super) type GeneName = FixedAscii<32>;

type GeneNames = Array1<GeneName>;

type FeatureType = FixedAscii<32>;

type FeatureTypes = Array1<FeatureType>;

#[derive(Debug, PartialEq)]
pub(crate) struct Features {
    ensembl_ids: EnsemblIds,
    gene_names: GeneNames,
    feature_types: FeatureTypes,
}

impl Features {
    pub(super) fn ensembl_ids(&self) -> &EnsemblIds {
        &self.ensembl_ids
    }

    pub(super) fn gene_names(&self) -> &GeneNames {
        &self.gene_names
    }

    pub(super) fn feature_types(&self) -> &FeatureTypes {
        &self.feature_types
    }

    pub(crate) fn len(&self) -> usize {
        self.ensembl_ids.len()
    }
}

#[derive(Clone, Serialize, Debug, thiserror::Error)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum VarError {
    #[error(
        "one or more fields in .var are missing or improperly formatted - ensure the correct \
         column names were provided"
    )]
    InvalidH5Fields { errors: Vec<ReadH5FieldError> },
    #[error("some genes were filtered out of the dataset (expected {}, found {n_found_genes}) - provide a dataset containing every gene in the transcriptome",
        n_expected_genes2.map_or_else(|| n_expected_genes.to_string(), |n2| format!("{n_expected_genes} or {n2}")))]
    FilteredGenes {
        n_expected_genes: usize,
        #[serde(skip_serializing_if = "Option::is_none")]
        n_expected_genes2: Option<usize>,
        n_found_genes: usize,
    },
    #[error(
        ".var has {ensembl_ids_len} Ensembl IDs, {gene_names_len} gene names and \
         {feature_types_len} feature types, but these must all be equal - regenerate the dataset \
         with scanpy"
    )]
    InvalidShapes {
        ensembl_ids_len: usize,
        gene_names_len: usize,
        feature_types_len: usize,
    },
    #[error(
        "{} genes in .var do not match the reference transcriptome - {}",
        errors.len(),
        hints(errors).join("; ")
    )]
    Genes { errors: Vec<VarRowError> },
}

#[derive(Debug, Clone, Serialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum VarRowError {
    DuplicateGene {
        ensembl_id: String,
        gene_name: String,
    },

    EnsemblIdGeneNameMismatch {
        ensembl_id: String,
        expected_gene_name: &'static str,
        found_gene_name: String,
    },
    UnrecognizedEnsemblId {
        ensembl_id: String,
    },
}

impl VarRowError {
    fn hint(&self) -> &'static str {
        match self {
            Self::DuplicateGene { .. } => "remove the duplicated genes from the dataset",
            Self::EnsemblIdGeneNameMismatch { .. } => {
                "if you used AnnData.var_names_make_unique, regenerate the dataset without it"
            }
            Self::UnrecognizedEnsemblId { .. } => {
                "ensure the dataset was aligned against the transcriptome you specified"
            }
        }
    }
}

#[cfg(test)]
mod tests {

    use hdf5_metno::{File, types::VarLenUnicode};
    use ndarray::arr1;

    use crate::{
        reference_dataset::{
            columns::{EnsemblIdCol, GeneNameCol},
            transcriptome::{Transcriptome, TranscriptomeName},
            var::{
                VarError, VarRowError, read_features_from_h5ad, validate_var_matches_transcriptome,
            },
        },
        target_list::chemistry::tests::tp53_ensembl_id,
    };

    fn read_generated_features(
        ensembl_id_col: &str,
        gene_name_col: &str,
    ) -> Result<super::Features, VarError> {
        read_features_from_h5ad(
            &File::open("test-data/csr_adata.h5ad").unwrap(),
            &EnsemblIdCol(ensembl_id_col.to_owned()),
            &GeneNameCol(gene_name_col.to_owned()),
            Transcriptome::new(TranscriptomeName::Grch382020A, false),
        )
    }

    #[test]
    fn filtered_genes_are_rejected() {
        // The generated datasets have 100 genes
        let error = read_generated_features("ensembl_id", "gene_name").unwrap_err();

        std::assert_matches!(
            error,
            VarError::FilteredGenes {
                n_found_genes: 100,
                ..
            }
        );
    }

    #[test]
    fn missing_var_columns_are_collected() {
        let VarError::InvalidH5Fields { errors } =
            read_generated_features("nonexistent", "also_nonexistent").unwrap_err()
        else {
            unreachable!();
        };

        // Only the two nonexistent columns should be reported - 'feature_types'
        // exists in the test-data
        assert_eq!(
            errors.len(),
            2,
            "every unreadable column should be reported, not just the first"
        );
    }

    #[test]
    fn all_gene_errors_are_collected() {
        let ensembl_ids = unsafe {
            [
                VarLenUnicode::from_str_unchecked(tp53_ensembl_id().as_str()),
                VarLenUnicode::from_str_unchecked("ENSG00000116678"),
            ]
        };

        let gene_names = unsafe {
            [
                VarLenUnicode::from_str_unchecked("foo"),
                VarLenUnicode::from_str_unchecked("bar"),
            ]
        };

        let transcriptome = Transcriptome::new(TranscriptomeName::Grch382024A, false).unwrap();

        let VarError::Genes { errors } = validate_var_matches_transcriptome(
            &arr1(&ensembl_ids),
            &arr1(&gene_names),
            transcriptome.gene_map(transcriptome.n_genes().0).unwrap(),
        )
        .unwrap_err() else {
            unreachable!();
        };

        std::assert_matches!(
            errors.as_slice(),
            [
                VarRowError::EnsemblIdGeneNameMismatch { .. },
                VarRowError::EnsemblIdGeneNameMismatch { .. }
            ]
        );
    }
}
