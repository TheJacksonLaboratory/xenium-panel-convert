use serde::Serialize;

use crate::target_list::{
    ValidTarget,
    target::{self, TargetId, TargetName, ValidGene},
};

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(transparent)]
pub struct XeniumPanelDesignerGeneList(Vec<XeniumPanelDesignerGene>);

impl XeniumPanelDesignerGeneList {
    pub fn from_valid_targets(mut valid_targets: Vec<ValidTarget>) -> Self {
        valid_targets.sort_by_key(super::target::ValidTarget::priority);

        Self(
            valid_targets
                .iter()
                .map(XeniumPanelDesignerGene::from_valid_target)
                .collect(),
        )
    }

    #[must_use]
    pub fn as_slice(&self) -> &[XeniumPanelDesignerGene] {
        &self.0
    }

    pub(crate) fn len(&self) -> usize {
        self.0.iter().len()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct XeniumPanelDesignerGene {
    #[serde(rename = "Gene")]
    gene: TargetName,
    #[serde(rename = "Ensembl ID")]
    ensembl_id: TargetId,
    #[serde(rename = "Probe sets")]
    probe_sets: Option<u16>,
    #[serde(rename = "Force")]
    force: Option<Force>,
}

impl XeniumPanelDesignerGene {
    fn from_valid_target(target: &ValidTarget) -> Self {
        let (ensembl_id, gene_name) = target.id_and_name();

        Self {
            gene: gene_name,
            ensembl_id,
            probe_sets: None,
            force: (target.priority() == target::Priority::MustHave).then_some(Force::Forced),
        }
    }

    pub(crate) fn gene(&self) -> Option<ValidGene> {
        let (TargetId::EnsemblId(ensembl_id), TargetName::GeneName(gene_name)) =
            (&self.ensembl_id, &self.gene)
        else {
            return None;
        };

        Some(ValidGene {
            ensembl_id: *ensembl_id,
            gene_name: *gene_name,
        })
    }
}

#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
enum Force {
    Forced,
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use crate::target_list::{
        chemistry::xenium_v1_human_ensembl_id_to_gene,
        parse_target_list,
        xenium_panel_designer::{Force, XeniumPanelDesignerGeneList},
    };

    fn gene_list() -> XeniumPanelDesignerGeneList {
        // Deliberately not in priority order
        let target_list = "ensembl_id,gene_name,group,priority\nENSG00000116678,LEPR,group0,\
                           backup\nENSG00000141510,TP53,group0,must_have\nENSG00000120802,TMPO,\
                           group1,desired";

        let targets = parse_target_list(
            target_list,
            &HashMap::new(),
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap();

        XeniumPanelDesignerGeneList::from_valid_targets(targets)
    }

    #[test]
    fn targets_are_sorted_by_priority() {
        let XeniumPanelDesignerGeneList(genes) = gene_list();

        let gene_names: Vec<_> = genes.iter().map(|g| g.gene.to_string()).collect();
        assert_eq!(
            gene_names,
            ["TP53", "TMPO", "LEPR"],
            "genes should be ordered must_have, desired, backup"
        );

        assert_eq!(genes[0].force, Some(Force::Forced));
        assert_eq!(
            genes[1].force, None,
            "only must_have targets should be forced"
        );
    }

    #[test]
    fn serializes_the_panel_designer_columns() {
        let mut writer = csv::Writer::from_writer(Vec::new());

        let XeniumPanelDesignerGeneList(genes) = gene_list();
        for gene in &genes {
            writer.serialize(gene).unwrap();
        }

        let csv = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        let mut rows = csv.lines();

        assert_eq!(rows.next().unwrap(), "Gene,Ensembl ID,Probe sets,Force");
    }
}
