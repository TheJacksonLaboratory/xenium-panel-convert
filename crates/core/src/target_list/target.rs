use std::str::FromStr;

use csv::StringRecord;
use serde::{Deserialize, Serialize};
use strum::VariantNames;

use crate::{
    common::ErrorVecExt,
    target_list::{
        TargetErrorInner,
        chemistry::{EnsemblId, GeneName, UnvalidatedEnsemblId, UnvalidatedGeneName},
    },
};

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct UnvalidatedGene {
    pub ensembl_id: Option<UnvalidatedEnsemblId>,
    pub gene_name: Option<UnvalidatedGeneName>,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct UnvalidatedTarget {
    #[serde(flatten)]
    pub gene: UnvalidatedGene,
    pub group: Option<String>,
    pub priority: Option<String>,
    pub custom: Option<String>,
}

impl UnvalidatedTarget {
    pub(super) fn from_record(record: &StringRecord, fieldnames: &StringRecord) -> Self {
        // Unwrapping is fine because extra fields won't cause a failure, nor will
        // missing fields
        record.deserialize(Some(fieldnames)).unwrap()
    }
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, Hash)]
pub struct ValidGene {
    pub(crate) ensembl_id: EnsemblId,
    pub(crate) gene_name: GeneName,
}

impl ValidGene {
    fn from_unvalidated(
        UnvalidatedGene {
            ensembl_id,
            gene_name: submitted_gene_name,
        }: &UnvalidatedGene,
        ensembl_id_to_gene: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)>,
    ) -> Result<Self, TargetErrorInner> {
        let Some(ensembl_id) = ensembl_id else {
            return Err(TargetErrorInner::NoEnsemblId);
        };

        let map_valid_gene = |(ensembl_id, gene_name)| Self {
            ensembl_id,
            gene_name,
        };

        if !ensembl_id.is_versionless_and_uppercase() {
            let correct_gene =
                ensembl_id_to_gene(&ensembl_id.to_versionless_uppercase()).map(map_valid_gene);

            return Err(TargetErrorInner::VersionedOrLowercaseEnsemblId { correct_gene });
        }

        let valid_gene = ensembl_id_to_gene(ensembl_id)
            .map(map_valid_gene)
            .ok_or(TargetErrorInner::GeneNotFound)?;

        let Some(submitted_gene_name) = submitted_gene_name else {
            return Err(TargetErrorInner::NoGeneName {
                probable_gene_name: valid_gene.gene_name,
            });
        };

        if *submitted_gene_name == valid_gene.gene_name {
            Ok(valid_gene)
        } else {
            Err(TargetErrorInner::EnsemblIdGeneNameMismatch {
                ensembl_id: valid_gene.ensembl_id,
                correct_gene_name: valid_gene.gene_name,
            })
        }
    }
}

#[derive(
    Clone,
    Copy,
    Debug,
    Serialize,
    PartialEq,
    Eq,
    strum::EnumString,
    PartialOrd,
    Ord,
    strum::VariantNames,
)]
#[serde(rename_all = "snake_case")]
#[strum(serialize_all = "snake_case")]
pub(super) enum Priority {
    MustHave,
    Desired,
    Backup,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq, Hash)]
pub(crate) enum TargetId {
    EnsemblId(EnsemblId),
    Custom(Option<UnvalidatedEnsemblId>),
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq, Hash)]
pub(crate) enum TargetName {
    GeneName(GeneName),
    Custom(Option<UnvalidatedGeneName>),
}

#[cfg(test)]
impl std::fmt::Display for TargetName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::GeneName(g) => g.fmt(f),
            Self::Custom(c) => c
                .as_ref()
                .map(UnvalidatedGeneName::as_str)
                .unwrap_or_default()
                .fmt(f),
        }
    }
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct ValidTarget {
    ensembl_id: TargetId,
    gene_name: TargetName,
    group: String,
    priority: Priority,
    custom: bool,
}

impl ValidTarget {
    pub(crate) fn id_and_name(&self) -> (TargetId, TargetName) {
        (self.ensembl_id.clone(), self.gene_name.clone())
    }

    pub(super) fn priority(&self) -> Priority {
        self.priority
    }

    pub(super) fn from_unvalidated(
        UnvalidatedTarget {
            gene,
            group,
            priority,
            custom,
        }: &UnvalidatedTarget,
        ensembl_id_to_gene: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)>,
    ) -> Result<Self, Vec<TargetErrorInner>> {
        let mut errors = Vec::new();

        let group = group.as_deref().map(str::to_ascii_lowercase);

        if group.is_none() {
            errors.push(TargetErrorInner::MissingField { fieldname: "group" });
        }

        let priority =
            parse_priority_field(priority.as_deref()).map_or_else(|err| errors.push_err(err), Some);

        let is_custom =
            parse_is_custom(custom.as_deref()).map_or_else(|err| errors.push_err(err), Some);

        let valid_gene = match is_custom {
            Some(true) => None,
            Some(false) => ValidGene::from_unvalidated(&gene, ensembl_id_to_gene)
                .map_or_else(|err| errors.push_err(err), Some),
            None => None,
        };

        match (valid_gene, group, priority, is_custom) {
            (
                Some(ValidGene {
                    ensembl_id,
                    gene_name,
                }),
                Some(group),
                Some(priority),
                Some(custom),
            ) => Ok(ValidTarget {
                ensembl_id: TargetId::EnsemblId(ensembl_id),
                gene_name: TargetName::GeneName(gene_name),
                group,
                priority,
                custom,
            }),
            (None, Some(group), Some(priority), Some(true)) => Ok(ValidTarget {
                ensembl_id: TargetId::Custom(gene.ensembl_id.clone()),
                gene_name: TargetName::Custom(gene.gene_name.clone()),
                group,
                priority,
                custom: true,
            }),
            _ => Err(errors),
        }
    }
}

fn parse_is_custom(s: Option<&str>) -> Result<bool, TargetErrorInner> {
    let Some(s) = s else {
        return Ok(false);
    };

    if s.eq_ignore_ascii_case("true") {
        Ok(true)
    } else if s.eq_ignore_ascii_case("false") {
        Ok(false)
    } else {
        Err(TargetErrorInner::MalformedCsv {
            reason: "field 'custom' must be 'true', 'false', or empty ('')".to_owned(),
        })
    }
}

fn parse_priority_field(s: Option<&str>) -> Result<Priority, TargetErrorInner> {
    let Some(s) = s else {
        return Err(TargetErrorInner::MissingField {
            fieldname: "priority",
        });
    };

    Priority::from_str(s).map_err(|_| TargetErrorInner::InvalidPriority {
        value: s.to_owned(),
        allowed: Priority::VARIANTS,
    })
}

#[cfg(test)]
mod tests {
    use strum::VariantNames;

    use crate::target_list::{
        TargetErrorInner,
        chemistry::{
            UnvalidatedEnsemblId, UnvalidatedGeneName, tests::tp53_ensembl_id,
            xenium_v1_human_ensembl_id_to_gene,
        },
        csv_util::read_csv_trimmed,
        target::{Priority, UnvalidatedGene, UnvalidatedTarget, ValidGene, ValidTarget},
    };

    #[test]
    fn valid_target() {
        let target = UnvalidatedTarget {
            gene: UnvalidatedGene {
                ensembl_id: Some(tp53_ensembl_id()),
                gene_name: Some(UnvalidatedGeneName::new("TP53".to_owned())),
            },
            group: Some("Group0".to_owned()),
            priority: Some("must_have".to_owned()),
            custom: None,
        };

        let valid_target =
            ValidTarget::from_unvalidated(&target, xenium_v1_human_ensembl_id_to_gene).unwrap();

        assert_eq!(valid_target.group, "group0", "group was not lowercased");
    }

    #[test]
    fn validate_collects_every_error_in_a_row() {
        let target = UnvalidatedTarget {
            gene: UnvalidatedGene {
                ensembl_id: None,
                gene_name: None,
            },
            group: None,
            priority: Some("urgent".to_owned()),
            custom: None,
        };

        let errors =
            ValidTarget::from_unvalidated(&target, xenium_v1_human_ensembl_id_to_gene).unwrap_err();

        assert_eq!(
            errors,
            [
                TargetErrorInner::MissingField { fieldname: "group" },
                TargetErrorInner::InvalidPriority {
                    value: "urgent".to_owned(),
                    allowed: Priority::VARIANTS,
                },
                TargetErrorInner::NoEnsemblId,
            ],
            "every error in a row should be reported, not just the first"
        );
    }

    #[test]
    fn no_ensembl_id() {
        let err = ValidGene::from_unvalidated(
            &UnvalidatedGene {
                ensembl_id: None,
                gene_name: Some(UnvalidatedGeneName::new("TP53".to_owned())),
            },
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        assert_eq!(err, TargetErrorInner::NoEnsemblId);
    }

    #[test]
    fn valid_gene() {
        ValidGene::from_unvalidated(
            &UnvalidatedGene {
                ensembl_id: Some(tp53_ensembl_id()),
                gene_name: Some(UnvalidatedGeneName::new("TP53".to_owned())),
            },
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap();
    }

    #[test]
    fn unvalidated_target_deserializes_with_invalid_fields() {
        let data = "field,ensembl_id\nvalue2,id";
        let mut reader = read_csv_trimmed(data);

        let fieldnames = reader.headers().unwrap().clone();
        let record = reader.records().next().unwrap().unwrap();
        let deserialized = UnvalidatedTarget::from_record(&record, &fieldnames);

        assert_eq!(
            deserialized,
            UnvalidatedTarget {
                gene: UnvalidatedGene {
                    ensembl_id: Some(UnvalidatedEnsemblId::new("id".to_owned())),
                    gene_name: None
                },
                group: None,
                priority: None,
                custom: None
            }
        );
    }

    #[test]
    fn ensembl_id_gene_name_mismatch() {
        let ensembl_id = tp53_ensembl_id();
        let gene_name = UnvalidatedGeneName::new(String::new());

        let err = ValidGene::from_unvalidated(
            &UnvalidatedGene {
                ensembl_id: Some(ensembl_id.clone()),
                gene_name: Some(gene_name.clone()),
            },
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        let (correct_ensembl_id, correct_gene_name) =
            xenium_v1_human_ensembl_id_to_gene(&ensembl_id).unwrap();

        assert_eq!(
            err,
            TargetErrorInner::EnsemblIdGeneNameMismatch {
                ensembl_id: correct_ensembl_id,
                correct_gene_name
            },
            "failed to create Ensembl ID-gene name mismatch error"
        );
    }

    #[test]
    fn versioned_or_lowercase_ensembl_id_suggests_correct_gene() {
        let ensembl_id = tp53_ensembl_id();
        let (correct_ensembl_id, correct_gene_name) =
            xenium_v1_human_ensembl_id_to_gene(&ensembl_id).unwrap();

        let versioned =
            UnvalidatedEnsemblId::new(format!("{}.1", ensembl_id.as_str().to_lowercase()));

        let err = ValidGene::from_unvalidated(
            &UnvalidatedGene {
                ensembl_id: Some(versioned),
                gene_name: Some(UnvalidatedGeneName::new("TP53".to_owned())),
            },
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        assert_eq!(
            err,
            TargetErrorInner::VersionedOrLowercaseEnsemblId {
                correct_gene: Some(ValidGene {
                    ensembl_id: correct_ensembl_id,
                    gene_name: correct_gene_name,
                }),
            }
        );
    }

    #[test]
    fn missing_gene_name_suggests_probable_name() {
        let ensembl_id = tp53_ensembl_id();
        let (_, correct_gene_name) = xenium_v1_human_ensembl_id_to_gene(&ensembl_id).unwrap();

        let err = ValidGene::from_unvalidated(
            &UnvalidatedGene {
                ensembl_id: Some(ensembl_id),
                gene_name: None,
            },
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        assert_eq!(
            err,
            TargetErrorInner::NoGeneName {
                probable_gene_name: correct_gene_name
            }
        );
    }
}
