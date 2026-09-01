#![allow(clippy::unreadable_literal)]
use std::fmt::Display;

use serde::{Deserialize, Serialize};

mod xenium_prime_human;
mod xenium_prime_mouse;
mod xenium_v1_human;
mod xenium_v1_mouse;

#[derive(Debug, Clone, Copy, strum::Display, serde::Serialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum Species {
    #[serde(rename = "Homo sapiens")]
    #[strum(serialize = "Homo sapiens")]
    HomoSapiens,
    #[serde(rename = "Mus musculus")]
    #[strum(serialize = "Mus musculus")]
    MusMusculus,
}

#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum Chemistry {
    V1,
    Prime,
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EnsemblId(&'static str);

impl EnsemblId {
    #[must_use]
    pub fn as_str(&self) -> &str {
        self.0
    }
}

impl Display for EnsemblId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GeneName(&'static str);

impl PartialEq<&str> for GeneName {
    fn eq(&self, other: &&str) -> bool {
        self.0 == *other
    }
}

impl Display for GeneName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct UnvalidatedEnsemblId(String);

impl UnvalidatedEnsemblId {
    #[cfg(test)]
    pub(super) fn new(ensembl_id: String) -> Self {
        Self(ensembl_id)
    }

    #[must_use]
    pub(super) fn to_versionless_uppercase(&self) -> Self {
        Self(self.0.split('.').next().unwrap_or("").to_uppercase())
    }

    #[must_use]
    pub(super) fn is_versionless_and_uppercase(&self) -> bool {
        !self
            .0
            .chars()
            .any(|c| c == '.' || (c.is_alphabetic() && !c.is_uppercase()))
    }

    #[cfg(test)]
    pub(crate) fn as_str(&self) -> &str {
        &self.0
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct UnvalidatedGeneName(String);

#[cfg(test)]
impl UnvalidatedGeneName {
    pub(super) fn new(gene_name: String) -> Self {
        Self(gene_name)
    }

    pub(super) fn as_str(&self) -> &str {
        &self.0
    }
}

impl PartialEq<GeneName> for UnvalidatedGeneName {
    fn eq(&self, other: &GeneName) -> bool {
        self.0 == other.0
    }
}

#[must_use]
pub fn xenium_v1_human_ensembl_id_to_gene(
    ensembl_id: &UnvalidatedEnsemblId,
) -> Option<(EnsemblId, GeneName)> {
    ensembl_id_to_gene(ensembl_id, &xenium_v1_human::XENIUM_V1_HUMAN_GENES)
}

#[must_use]
pub fn xenium_prime_human_ensembl_id_to_gene(
    ensembl_id: &UnvalidatedEnsemblId,
) -> Option<(EnsemblId, GeneName)> {
    ensembl_id_to_gene(ensembl_id, &xenium_prime_human::XENIUM_PRIME_HUMAN_GENES)
}

#[must_use]
pub fn xenium_v1_mouse_ensembl_id_to_gene(
    ensembl_id: &UnvalidatedEnsemblId,
) -> Option<(EnsemblId, GeneName)> {
    ensembl_id_to_gene(ensembl_id, &xenium_v1_mouse::XENIUM_V1_MOUSE_GENES)
}

#[must_use]
pub fn xenium_prime_mouse_ensembl_id_to_gene(
    ensembl_id: &UnvalidatedEnsemblId,
) -> Option<(EnsemblId, GeneName)> {
    ensembl_id_to_gene(ensembl_id, &xenium_prime_mouse::XENIUM_PRIME_MOUSE_GENES)
}

fn ensembl_id_to_gene(
    ensembl_id: &UnvalidatedEnsemblId,
    map: &phf::Map<&'static str, &'static str>,
) -> Option<(EnsemblId, GeneName)> {
    map.get_entry(&ensembl_id.0)
        .map(|(eid, gn)| (EnsemblId(eid), GeneName(gn)))
}

#[cfg(test)]
pub(crate) mod tests {
    use crate::target_list::chemistry::{
        GeneName, UnvalidatedEnsemblId, xenium_prime_human::XENIUM_PRIME_HUMAN_GENES,
        xenium_prime_mouse::XENIUM_PRIME_MOUSE_GENES, xenium_v1_human::XENIUM_V1_HUMAN_GENES,
        xenium_v1_human_ensembl_id_to_gene, xenium_v1_mouse::XENIUM_V1_MOUSE_GENES,
    };

    pub(crate) fn tp53_ensembl_id() -> UnvalidatedEnsemblId {
        UnvalidatedEnsemblId("ENSG00000141510".to_owned())
    }

    #[test]
    fn unavailable_genes_are_not_in_map() {
        let unavailable_genes = [
            ("ENSG00000273816", &XENIUM_V1_HUMAN_GENES),
            ("ENSG00000249966", &XENIUM_PRIME_HUMAN_GENES),
            ("ENSMUSG00000117061", &XENIUM_V1_MOUSE_GENES),
            ("ENSMUSG00000094028", &XENIUM_PRIME_MOUSE_GENES),
        ];

        for (ensembl_id, map) in unavailable_genes {
            std::assert_matches!(map.get(ensembl_id), None);
        }
    }

    #[test]
    fn correct_detection_of_versionless_uppercase_ensembl_id() {
        let ensembl_id = tp53_ensembl_id();
        assert!(ensembl_id.is_versionless_and_uppercase());

        let ensembl_id = UnvalidatedEnsemblId(format!("{}.1", ensembl_id.0.to_lowercase()));
        assert!(!ensembl_id.is_versionless_and_uppercase());
    }

    #[test]
    fn canonicalized_ensembl_id_gets_correct_gene_name() {
        let ensembl_id = tp53_ensembl_id().to_versionless_uppercase();

        std::assert_matches!(
            xenium_v1_human_ensembl_id_to_gene(&ensembl_id).unwrap(),
            (_, GeneName("TP53"))
        );
    }
}
