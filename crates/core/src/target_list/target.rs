use std::{collections::HashMap, str::FromStr};

use csv::StringRecord;
use serde::{Deserialize, Serialize};
use strum::VariantNames;

use crate::{
    error::collect_error,
    target_list::{
        TargetError, TargetList,
        chemistry::{EnsemblId, GeneSymbol, UnvalidatedEnsemblId, UnvalidatedGeneName},
    },
};

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq, Hash)]
pub struct UnvalidatedGene {
    pub ensembl_id: Option<UnvalidatedEnsemblId>,
    pub gene_symbol: Option<UnvalidatedGeneName>,
}

pub(super) const FIELDNAMES: [&str; 5] = ["ensembl_id", "gene_symbol", "group", "priority", "custom"];

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct UnvalidatedTarget {
    #[serde(flatten)]
    pub gene: UnvalidatedGene,
    pub group: Option<String>,
    pub priority: Option<String>,
    pub custom: Option<String>,
    #[serde(flatten)]
    pub other_fields: HashMap<String, serde_json::Value>,
}

impl UnvalidatedTarget {
    pub(super) fn from_record(record: &StringRecord, fieldnames: &StringRecord) -> Self {
        // Unwrapping is fine because extra fields won't cause a failure, nor
        // will missing fields
        record.deserialize(Some(fieldnames)).unwrap()
    }
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, Hash)]
pub struct ValidGene {
    pub(crate) ensembl_id: EnsemblId,
    pub(crate) gene_symbol: GeneSymbol,
}

impl ValidGene {
    fn from_unvalidated(
        UnvalidatedGene {
            ensembl_id,
            gene_symbol: submitted_gene_symbol,
        }: &UnvalidatedGene,
        ensembl_id_to_gene: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneSymbol)>,
    ) -> Result<Self, TargetError> {
        let Some(ensembl_id) = ensembl_id else {
            return Err(TargetError::NoEnsemblId);
        };

        let map_valid_gene = |(ensembl_id, gene_symbol)| Self {
            ensembl_id,
            gene_symbol,
        };

        if !ensembl_id.is_versionless_and_uppercase() {
            let correct_gene = ensembl_id_to_gene(&ensembl_id.to_versionless_uppercase()).map(map_valid_gene);

            return Err(TargetError::VersionedOrLowercaseEnsemblId { correct_gene });
        }

        let valid_gene = ensembl_id_to_gene(ensembl_id)
            .map(map_valid_gene)
            .ok_or(TargetError::GeneNotFound)?;

        let Some(submitted_gene_symbol) = submitted_gene_symbol else {
            return Err(TargetError::NoGeneName {
                probable_gene_symbol: valid_gene.gene_symbol,
            });
        };

        if *submitted_gene_symbol == valid_gene.gene_symbol {
            Ok(valid_gene)
        } else {
            Err(TargetError::EnsemblIdGeneNameMismatch {
                ensembl_id: valid_gene.ensembl_id,
                correct_gene_symbol: valid_gene.gene_symbol,
            })
        }
    }
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, strum::EnumString, PartialOrd, Ord, strum::VariantNames)]
#[serde(rename_all = "snake_case")]
#[strum(serialize_all = "snake_case")]
pub(super) enum Priority {
    MustHave,
    Desired,
    Backup,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize)]
pub(crate) enum TargetGene {
    Standard(ValidGene),
    Custom(UnvalidatedGene),
}

impl TargetGene {
    pub(crate) fn as_strs(&self) -> (Option<&str>, Option<&str>) {
        match self {
            Self::Standard(gene) => (Some(gene.ensembl_id.as_str()), Some(gene.gene_symbol.as_str())),
            Self::Custom(gene) => (
                gene.ensembl_id.as_ref().map(UnvalidatedEnsemblId::as_str),
                gene.gene_symbol.as_ref().map(UnvalidatedGeneName::as_str),
            ),
        }
    }

    pub(crate) fn valid_gene(&self) -> Option<ValidGene> {
        match self {
            Self::Standard(gene) => Some(*gene),
            Self::Custom(_) => None,
        }
    }

    fn is_custom(&self) -> bool {
        matches!(self, Self::Custom(_))
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ValidTarget {
    gene: TargetGene,
    group: String,
    priority: Priority,
    other_fields: HashMap<String, serde_json::Value>,
}

impl ValidTarget {
    pub(crate) fn gene(&self) -> &TargetGene {
        &self.gene
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
            other_fields,
        }: &UnvalidatedTarget,
        ensembl_id_to_gene: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneSymbol)>,
    ) -> Result<Self, Vec<TargetError>> {
        let mut errors = Vec::new();

        let group = group.as_deref().map(str::to_ascii_lowercase);

        if group.is_none() {
            errors.push(TargetError::MissingField { fieldname: "group" });
        }

        let priority = collect_error(parse_priority_field(priority.as_deref()), &mut errors);

        let is_custom = collect_error(parse_custom_field(custom.as_deref()), &mut errors);

        let target_gene = match is_custom {
            Some(false) => collect_error(ValidGene::from_unvalidated(gene, ensembl_id_to_gene), &mut errors)
                .map(TargetGene::Standard),
            Some(true) => Some(TargetGene::Custom(gene.clone())),
            None => None,
        };

        match (target_gene, group, priority) {
            (Some(gene), Some(group), Some(priority)) => Ok(ValidTarget {
                gene,
                group,
                priority,
                other_fields: other_fields.clone(),
            }),
            _ => Err(errors),
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct ValidTargetCsvRow<'a> {
    ensembl_id: Option<&'a str>,
    gene_symbol: Option<&'a str>,
    group: &'a str,
    priority: Priority,
    custom: bool,
    other_values: Vec<serde_json::Value>,
}

impl<'a> ValidTargetCsvRow<'a> {
    fn from_valid_target(
        ValidTarget {
            gene,
            group,
            priority,
            other_fields,
        }: &'a ValidTarget,
        custom_fieldnames: &[String],
    ) -> Self {
        let (ensembl_id, gene_symbol) = gene.as_strs();

        let other_values = custom_fieldnames
            .iter()
            .map(|fieldname| {
                // Unwrapping is fine because every record has every column of
                // the header, and every column that isn't one
                // of FIELDNAMES is captured in `other_fields`
                other_fields
                    .get(fieldname)
                    .cloned()
                    .expect("custom field is missing from target")
            })
            .collect();

        Self {
            ensembl_id,
            gene_symbol,
            group,
            priority: *priority,
            custom: gene.is_custom(),
            other_values,
        }
    }
}

#[derive(Debug)]
pub struct ValidTargetCsv<'a> {
    header: Vec<&'a str>,
    rows: Vec<ValidTargetCsvRow<'a>>,
}

impl<'a> ValidTargetCsv<'a> {
    pub fn from_target_list(
        TargetList {
            targets,
            custom_fieldnames,
        }: &'a TargetList,
    ) -> Self {
        let header = FIELDNAMES
            .into_iter()
            .chain(custom_fieldnames.iter().map(String::as_str))
            .collect();

        let rows = targets
            .iter()
            .map(|target| ValidTargetCsvRow::from_valid_target(target, custom_fieldnames))
            .collect();

        Self { header, rows }
    }

    #[must_use]
    pub fn header(&self) -> &[&'a str] {
        &self.header
    }

    #[must_use]
    pub fn rows(&self) -> &[ValidTargetCsvRow<'a>] {
        &self.rows
    }
}

fn parse_custom_field(s: Option<&str>) -> Result<bool, TargetError> {
    let Some(s) = s else {
        return Ok(false);
    };

    if s.eq_ignore_ascii_case("true") {
        Ok(true)
    } else if s.eq_ignore_ascii_case("false") {
        Ok(false)
    } else {
        Err(TargetError::InvalidValue {
            field: "custom",
            value: s.to_owned(),
            allowed: &["true", "false"],
        })
    }
}

fn parse_priority_field(s: Option<&str>) -> Result<Priority, TargetError> {
    let Some(s) = s else {
        return Err(TargetError::MissingField { fieldname: "priority" });
    };

    Priority::from_str(s).map_err(|_| TargetError::InvalidValue {
        field: "priority",
        value: s.to_owned(),
        allowed: Priority::VARIANTS,
    })
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use strum::VariantNames;

    use crate::target_list::{
        TargetError,
        chemistry::{
            UnvalidatedEnsemblId, UnvalidatedGeneName, tests::tp53_ensembl_id, xenium_v1_human_ensembl_id_to_gene,
        },
        csv_util::read_csv_trimmed,
        parse_target_list,
        target::{
            FIELDNAMES, Priority, UnvalidatedGene, UnvalidatedTarget, ValidGene, ValidTarget, ValidTargetCsv,
            ValidTargetCsvRow,
        },
    };

    fn serialize_csv(ValidTargetCsv { header, rows }: &ValidTargetCsv<'_>) -> Vec<u8> {
        let mut writer = csv::WriterBuilder::new().has_headers(false).from_writer(Vec::new());

        writer.write_record(header).unwrap();

        for row in rows {
            writer.serialize(row).unwrap();
        }

        writer.into_inner().unwrap()
    }

    #[test]
    fn valid_target() {
        let target = UnvalidatedTarget {
            gene: UnvalidatedGene {
                ensembl_id: Some(tp53_ensembl_id()),
                gene_symbol: Some(UnvalidatedGeneName::new("TP53".to_owned())),
            },
            group: Some("Group0".to_owned()),
            priority: Some("must_have".to_owned()),
            custom: None,
            other_fields: HashMap::new(),
        };

        let valid_target = ValidTarget::from_unvalidated(&target, xenium_v1_human_ensembl_id_to_gene).unwrap();

        assert_eq!(valid_target.group, "group0", "group was not lowercased");
    }

    #[test]
    fn validate_collects_every_error_in_a_row() {
        let target = UnvalidatedTarget {
            gene: UnvalidatedGene {
                ensembl_id: None,
                gene_symbol: None,
            },
            group: None,
            priority: Some("urgent".to_owned()),
            custom: None,
            other_fields: HashMap::new(),
        };

        let errors = ValidTarget::from_unvalidated(&target, xenium_v1_human_ensembl_id_to_gene).unwrap_err();

        assert_eq!(
            errors,
            [
                TargetError::MissingField { fieldname: "group" },
                TargetError::InvalidValue {
                    field: "priority",
                    value: "urgent".to_owned(),
                    allowed: Priority::VARIANTS,
                },
                TargetError::NoEnsemblId,
            ],
            "every error in a row should be reported, not just the first"
        );
    }

    #[test]
    fn no_ensembl_id() {
        let err = ValidGene::from_unvalidated(
            &UnvalidatedGene {
                ensembl_id: None,
                gene_symbol: Some(UnvalidatedGeneName::new("TP53".to_owned())),
            },
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        assert_eq!(err, TargetError::NoEnsemblId);
    }

    #[test]
    fn valid_gene() {
        ValidGene::from_unvalidated(
            &UnvalidatedGene {
                ensembl_id: Some(tp53_ensembl_id()),
                gene_symbol: Some(UnvalidatedGeneName::new("TP53".to_owned())),
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
                    gene_symbol: None
                },
                group: None,
                priority: None,
                custom: None,
                other_fields: [("field".to_owned(), "value2".into())].into_iter().collect(),
            }
        );
    }

    #[test]
    fn ensembl_id_gene_symbol_mismatch() {
        let ensembl_id = tp53_ensembl_id();
        let gene_symbol = UnvalidatedGeneName::new(String::new());

        let err = ValidGene::from_unvalidated(
            &UnvalidatedGene {
                ensembl_id: Some(ensembl_id.clone()),
                gene_symbol: Some(gene_symbol.clone()),
            },
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        let (correct_ensembl_id, correct_gene_symbol) = xenium_v1_human_ensembl_id_to_gene(&ensembl_id).unwrap();

        assert_eq!(
            err,
            TargetError::EnsemblIdGeneNameMismatch {
                ensembl_id: correct_ensembl_id,
                correct_gene_symbol
            },
            "failed to create Ensembl ID-gene name mismatch error"
        );
    }

    #[test]
    fn versioned_or_lowercase_ensembl_id_suggests_correct_gene() {
        let ensembl_id = tp53_ensembl_id();
        let (correct_ensembl_id, correct_gene_symbol) = xenium_v1_human_ensembl_id_to_gene(&ensembl_id).unwrap();

        let versioned = UnvalidatedEnsemblId::new(format!("{}.1", ensembl_id.as_str().to_lowercase()));

        let err = ValidGene::from_unvalidated(
            &UnvalidatedGene {
                ensembl_id: Some(versioned),
                gene_symbol: Some(UnvalidatedGeneName::new("TP53".to_owned())),
            },
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        assert_eq!(
            err,
            TargetError::VersionedOrLowercaseEnsemblId {
                correct_gene: Some(ValidGene {
                    ensembl_id: correct_ensembl_id,
                    gene_symbol: correct_gene_symbol,
                }),
            }
        );
    }

    #[test]
    fn missing_gene_symbol_suggests_probable_name() {
        let ensembl_id = tp53_ensembl_id();
        let (_, correct_gene_symbol) = xenium_v1_human_ensembl_id_to_gene(&ensembl_id).unwrap();

        let err = ValidGene::from_unvalidated(
            &UnvalidatedGene {
                ensembl_id: Some(ensembl_id),
                gene_symbol: None,
            },
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        assert_eq!(
            err,
            TargetError::NoGeneName {
                probable_gene_symbol: correct_gene_symbol
            }
        );
    }

    #[test]
    fn valid_target_serializes() {
        let row = ValidTargetCsvRow {
            ensembl_id: Some("some_ensembl_id"),
            gene_symbol: Some("some_gene_symbol"),
            group: "some_group",
            priority: Priority::MustHave,
            custom: true,
            other_values: vec![serde_json::Value::from("value")],
        };

        let csv = ValidTargetCsv {
            header: FIELDNAMES.into_iter().chain(["field"]).collect(),
            rows: vec![row],
        };

        let data = serialize_csv(&csv);

        assert_eq!(data, b"ensembl_id,gene_symbol,group,priority,custom,field\nsome_ensembl_id,some_gene_symbol,some_group,must_have,true,value\n");
    }

    #[test]
    fn custom_fields_are_propagated_in_input_order_with_correct_types() {
        let ensembl_id = tp53_ensembl_id();
        let ensembl_id = ensembl_id.as_str();
        let target_list = format!(
            "ensembl_id,str_field,gene_symbol,group,priority,number_field\n{ensembl_id},str_value,TP53,group0,\
             must_have,0"
        );

        let target_list = parse_target_list(&target_list, &HashMap::new(), xenium_v1_human_ensembl_id_to_gene).unwrap();

        assert_eq!(
            target_list.targets[0].other_fields["number_field"],
            serde_json::Value::from(0)
        );

        let data = serialize_csv(&ValidTargetCsv::from_target_list(&target_list));

        let expected = format!(
            "ensembl_id,gene_symbol,group,priority,custom,str_field,number_field\n{ensembl_id},TP53,group0,must_have,\
             false,str_value,0\n"
        );
        assert_eq!(data, expected.as_bytes());
    }
}
