use serde::Serialize;

use crate::target_list::target::{Priority, ValidTarget};

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct XeniumPanelDesignerCsv<'a> {
    header: [&'static str; 4],
    rows: Vec<XeniumPanelDesignerGene<'a>>,
}

impl<'a> XeniumPanelDesignerCsv<'a> {
    #[must_use]
    pub fn from_valid_targets(targets: &'a [ValidTarget]) -> Self {
        let mut targets_by_priority: Vec<&ValidTarget> = targets.iter().collect();
        targets_by_priority.sort_by_key(|target| target.priority());

        let rows = targets_by_priority
            .into_iter()
            .map(XeniumPanelDesignerGene::from_valid_target)
            .collect();

        Self {
            header: ["Gene", "Ensembl ID", "Probe sets", "Force"],
            rows,
        }
    }

    #[must_use]
    pub fn header(&self) -> &[&'static str] {
        &self.header
    }

    #[must_use]
    pub fn rows(&self) -> &[XeniumPanelDesignerGene<'a>] {
        &self.rows
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct XeniumPanelDesignerGene<'a> {
    #[serde(rename = "Gene")]
    gene: Option<&'a str>,
    #[serde(rename = "Ensembl ID")]
    ensembl_id: Option<&'a str>,
    #[serde(rename = "Probe sets")]
    probe_sets: Option<u16>,
    #[serde(rename = "Force")]
    force: Option<Force>,
}

impl<'a> XeniumPanelDesignerGene<'a> {
    fn from_valid_target(target: &'a ValidTarget) -> Self {
        let (ensembl_id, gene) = target.gene().as_strs();

        Self {
            gene,
            ensembl_id,
            probe_sets: None,
            force: (target.priority() == Priority::MustHave).then_some(Force::Forced),
        }
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
        target::ValidTarget,
        xenium_panel_designer::{Force, XeniumPanelDesignerCsv},
    };

    fn valid_targets() -> Vec<ValidTarget> {
        // Deliberately not in priority order
        let target_list = "ensembl_id,gene_symbol,group,priority\nENSG00000116678,LEPR,group0,\
                           backup\nENSG00000141510,TP53,group0,must_have\nENSG00000120802,TMPO,\
                           group1,desired";

        parse_target_list(
            target_list,
            &HashMap::new(),
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap()
        .targets
    }

    #[test]
    fn targets_are_sorted_by_priority() {
        let targets = valid_targets();
        let xpd_csv = XeniumPanelDesignerCsv::from_valid_targets(&targets);
        let genes = xpd_csv.rows();

        let gene_symbols: Vec<_> = genes.iter().map(|g| g.gene.unwrap()).collect();
        assert_eq!(
            gene_symbols,
            ["TP53", "TMPO", "LEPR"],
            "genes should be ordered must_have, desired, backup"
        );

        assert_eq!(genes[0].force, Some(Force::Forced));
        assert_eq!(
            genes[1].force, None,
            "only must_have targets should be forced"
        );
    }
}
