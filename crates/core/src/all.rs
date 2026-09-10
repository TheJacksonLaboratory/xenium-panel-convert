use crate::{
    error::Hinted,
    reference_dataset::{
        pseudo_anndata::PseudoAnndata,
        transcriptome::{Transcriptome, TranscriptomeName},
    },
    target_list::{
        chemistry::Species,
        target::{ValidGene, ValidTarget},
    },
};

#[must_use]
pub fn validate_target_list_and_reference_dataset_compatibility(
    target_list: &[ValidTarget],
    target_list_species: Species,
    reference_dataset: &PseudoAnndata,
    reference_dataset_transcriptome: TranscriptomeName,
    reference_dataset_is_flex: bool,
) -> Vec<Hinted<TargetListReferenceDatasetCompatibilityWarning>> {
    // Return early because this warning implies the other two warning types
    if !species_and_transcriptome_match(target_list_species, reference_dataset_transcriptome) {
        return vec![Hinted::new(
            TargetListReferenceDatasetCompatibilityWarning::SpeciesTranscriptomeMismatch {
                target_list_species,
                reference_dataset_transcriptome,
            },
        )];
    }

    let mut warnings = Vec::with_capacity(target_list.len());

    for target in target_list {
        let Some(gene) = target.gene().valid_gene() else {
            continue;
        };

        match validate_gene_is_in_transcriptome_with_correct_name(
            gene,
            reference_dataset,
            reference_dataset_transcriptome,
            reference_dataset_is_flex,
        ) {
            Ok(()) => (),
            Err(w) => {
                warnings.push(Hinted::new(w));
            }
        }
    }

    warnings
}

fn species_and_transcriptome_match(species: Species, transcriptome: TranscriptomeName) -> bool {
    matches!(
        (species, transcriptome),
        (
            Species::HomoSapiens,
            TranscriptomeName::Grch382020A | TranscriptomeName::Grch382024A
        ) | (
            Species::MusMusculus,
            TranscriptomeName::Mm102020A | TranscriptomeName::Grcm392024A
        )
    )
}

fn validate_gene_is_in_transcriptome_with_correct_name(
    gene: ValidGene,
    reference_dataset: &PseudoAnndata,
    reference_dataset_transcriptome: TranscriptomeName,
    reference_dataset_is_flex: bool,
) -> Result<(), TargetListReferenceDatasetCompatibilityWarning> {
    // If we have a PseudoAnndata, we know that the either the transcriptome is
    // 'other' or the features match the transcriptome exactly, so it's okay to
    // return Ok with no transcriptome
    let Some(transcriptome) = Transcriptome::new(reference_dataset_transcriptome, reference_dataset_is_flex) else {
        return Ok(());
    };

    let gene_map = transcriptome
        .gene_map(reference_dataset.features().len())
        .expect("if we have a PseudoAnndata, we know its features are exactly the transcriptome");

    let gene_symbol_from_transcriptome = gene_map.get(gene.ensembl_id.as_str()).ok_or(
        TargetListReferenceDatasetCompatibilityWarning::TargetNotInReferenceDataset {
            gene,
            transcriptome: reference_dataset_transcriptome,
        },
    )?;

    if gene.gene_symbol != *gene_symbol_from_transcriptome {
        return Err(TargetListReferenceDatasetCompatibilityWarning::GeneNameMismatch {
            gene_in_target_list: gene,
            gene_symbol_in_reference_dataset: gene_symbol_from_transcriptome,
        });
    }

    Ok(())
}

#[derive(Clone, Copy, Debug, thiserror::Error, serde::Serialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum TargetListReferenceDatasetCompatibilityWarning {
    #[error(
        "the target-list is {target_list_species} but the reference dataset was aligned against \
         {reference_dataset_transcriptome}"
    )]
    SpeciesTranscriptomeMismatch {
        target_list_species: Species,
        reference_dataset_transcriptome: TranscriptomeName,
    },
    #[error("{} ({}) is not in the reference dataset, whose transcriptome is {transcriptome}", gene.gene_symbol, gene.ensembl_id)]
    TargetNotInReferenceDataset {
        gene: ValidGene,
        transcriptome: TranscriptomeName,
    },
    #[error("{} is called {} in the target-list, but it is called {gene_symbol_in_reference_dataset} in the reference dataset - this is likely because the reference-dataset was aligned against GRCh38-2024-A or GRCm39-2024-A - add a new column to .var in the reference dataset with gene names from the 2020-A version of the transcriptome", gene_in_target_list.ensembl_id, gene_in_target_list.gene_symbol)]
    GeneNameMismatch {
        gene_in_target_list: ValidGene,
        gene_symbol_in_reference_dataset: &'static str,
    },
}
