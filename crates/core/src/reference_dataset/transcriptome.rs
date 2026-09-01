#![allow(clippy::unreadable_literal)]
use crate::reference_dataset::transcriptome::{
    grch38_2020_a::GRCH38_2020_A, grch38_2020_a_flex::GRCH38_2020_A_FLEX,
    grch38_2024_a::GRCH38_2024_A, grch38_2024_a_flex_v1_1::GRCH38_2024_A_FLEX_V1_1,
    grch38_2024_a_flex_v2_0::GRCH38_2024_A_FLEX_V2_0, grcm39_2024_a::GRCM39_2024_A,
    grcm39_2024_a_flex_v1_1::GRCM39_2024_A_FLEX_V1_1,
    grcm39_2024_a_flex_v2_0::GRCM39_2024_A_FLEX_V2_0, mm10_2020_a::MM10_2020_A,
    mm10_2020_a_flex::MM10_2020_A_FLEX,
};

mod grch38_2020_a;
mod grch38_2020_a_flex;
mod grch38_2024_a;
mod grch38_2024_a_flex_v1_1;
mod grch38_2024_a_flex_v2_0;
mod grcm39_2024_a;
mod grcm39_2024_a_flex_v1_1;
mod grcm39_2024_a_flex_v2_0;
mod mm10_2020_a;
mod mm10_2020_a_flex;

#[derive(Debug, Clone, Copy)]
pub enum Transcriptome {
    ThreePrime(&'static phf::Map<&'static str, &'static str>),
    Flex2020A(&'static phf::Map<&'static str, &'static str>),
    Flex2024A {
        v1: &'static phf::Map<&'static str, &'static str>,
        v2: &'static phf::Map<&'static str, &'static str>,
    },
}

impl Transcriptome {
    #[must_use]
    pub fn new(transcriptome: TranscriptomeName, flex: bool) -> Option<Self> {
        let transcriptome = match (transcriptome, flex) {
            (TranscriptomeName::Other, _) => {
                return None;
            }
            (TranscriptomeName::Grch382020A, false) => Self::ThreePrime(&GRCH38_2020_A),
            (TranscriptomeName::Grch382020A, true) => Self::Flex2020A(&GRCH38_2020_A_FLEX),
            (TranscriptomeName::Grch382024A, false) => Self::ThreePrime(&GRCH38_2024_A),
            (TranscriptomeName::Grch382024A, true) => Self::Flex2024A {
                v1: &GRCH38_2024_A_FLEX_V1_1,
                v2: &GRCH38_2024_A_FLEX_V2_0,
            },
            (TranscriptomeName::Mm102020A, false) => Self::ThreePrime(&MM10_2020_A),
            (TranscriptomeName::Mm102020A, true) => Self::Flex2020A(&MM10_2020_A_FLEX),
            (TranscriptomeName::Grcm392024A, false) => Self::ThreePrime(&GRCM39_2024_A),
            (TranscriptomeName::Grcm392024A, true) => Self::Flex2024A {
                v1: &GRCM39_2024_A_FLEX_V1_1,
                v2: &GRCM39_2024_A_FLEX_V2_0,
            },
        };

        Some(transcriptome)
    }

    #[must_use]
    pub(crate) fn n_genes(&self) -> (usize, Option<usize>) {
        match self {
            Self::ThreePrime(genes) | Self::Flex2020A(genes) => (genes.len(), None),
            Self::Flex2024A { v1, v2 } if v1.len() == v2.len() => (v1.len(), None),
            Self::Flex2024A { v1, v2 } => (v1.len(), Some(v2.len())),
        }
    }

    #[must_use]
    pub(crate) fn gene_map(
        self,
        n_genes_in_dataset: usize,
    ) -> Option<&'static phf::Map<&'static str, &'static str>> {
        match self {
            Self::ThreePrime(genes) | Self::Flex2020A(genes)
                if n_genes_in_dataset == genes.len() =>
            {
                Some(genes)
            }
            Self::Flex2024A { v1, v2 } => {
                let is_v1 = v1.len() == n_genes_in_dataset;
                let is_v2 = v2.len() == n_genes_in_dataset;

                let v1 = is_v1.then_some(v1);
                let v2 = is_v2.then_some(v2);

                // The only transcriptomes with the same number of genes are probe-sets
                // GRCm39-2024-A v1.1.1 and probe-sets GRCm39-2024-A v2.0.0. Luckily, they are
                // exactly the same, so just default to v2 since that's more recent
                v2.or(v1)
            }
            Self::ThreePrime(_) | Self::Flex2020A(_) => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, strum::EnumString, strum::Display, serde::Serialize)]
#[strum(ascii_case_insensitive)]
pub enum TranscriptomeName {
    #[serde(rename = "GRCh38-2020-A")]
    #[strum(serialize = "GRCh38-2020-A", serialize = "h2020")]
    Grch382020A,
    #[serde(rename = "GRCh38-2024-A")]
    #[strum(serialize = "GRCh38-2024-A", serialize = "h2024")]
    Grch382024A,
    #[serde(rename = "mm10-2020-A")]
    #[strum(serialize = "mm10-2020-A", serialize = "m2020")]
    Mm102020A,
    #[serde(rename = "GRCm39-2024-A")]
    #[strum(serialize = "GRCm39-2024-A", serialize = "m2024")]
    Grcm392024A,
    #[strum(serialize = "other", serialize = "o")]
    Other,
}

#[cfg(test)]
mod tests {
    use crate::reference_dataset::transcriptome::{
        Transcriptome, TranscriptomeName, grch38_2020_a::GRCH38_2020_A,
        grch38_2024_a_flex_v1_1::GRCH38_2024_A_FLEX_V1_1,
        grch38_2024_a_flex_v2_0::GRCH38_2024_A_FLEX_V2_0,
        grcm39_2024_a_flex_v2_0::GRCM39_2024_A_FLEX_V2_0,
    };

    #[test]
    fn genes_are_selected_by_count() {
        let three_prime = Transcriptome::new(TranscriptomeName::Grch382020A, false).unwrap();

        std::assert_matches!(three_prime.gene_map(100), None);
        assert_eq!(
            three_prime.gene_map(GRCH38_2020_A.len()).unwrap(),
            &GRCH38_2020_A
        );

        let human_flex = Transcriptome::new(TranscriptomeName::Grch382024A, true).unwrap();

        assert_eq!(
            human_flex.gene_map(GRCH38_2024_A_FLEX_V1_1.len()).unwrap(),
            &GRCH38_2024_A_FLEX_V1_1
        );
        assert_eq!(
            human_flex.gene_map(GRCH38_2024_A_FLEX_V2_0.len()).unwrap(),
            &GRCH38_2024_A_FLEX_V2_0
        );

        let mouse_flex = Transcriptome::new(TranscriptomeName::Grcm392024A, true).unwrap();

        assert_eq!(
            mouse_flex.gene_map(GRCM39_2024_A_FLEX_V2_0.len()).unwrap(),
            &GRCM39_2024_A_FLEX_V2_0,
            "expected v2.0.0 when both probe-sets have the same number of genes"
        );
    }
}
