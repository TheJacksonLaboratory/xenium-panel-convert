use std::fmt::Display;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CellBarcodeCol(pub String);

impl Default for CellBarcodeCol {
    fn default() -> Self {
        Self(String::from("_index"))
    }
}

impl Display for CellBarcodeCol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CellAnnotationCol(pub String);

impl Display for CellAnnotationCol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EnsemblIdCol(pub String);

impl EnsemblIdCol {
    #[must_use]
    pub(super) fn as_str(&self) -> &str {
        &self.0
    }
}

impl Display for EnsemblIdCol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl Default for EnsemblIdCol {
    fn default() -> Self {
        Self(String::from("gene_ids"))
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeneNameCol(pub String);

impl GeneNameCol {
    #[must_use]
    pub(super) fn as_str(&self) -> &str {
        &self.0
    }
}

impl Display for GeneNameCol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl Default for GeneNameCol {
    fn default() -> Self {
        Self(String::from("_index"))
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CountsLayerName(String);

impl CountsLayerName {
    #[must_use]
    pub(super) fn as_str(&self) -> &str {
        &self.0
    }

    pub fn new(layer_name: &str) -> Self {
        Self(format!("layers/{layer_name}"))
    }

    pub fn x() -> Self {
        Self::default()
    }
}

impl Display for CountsLayerName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl Default for CountsLayerName {
    fn default() -> Self {
        Self(String::from("X"))
    }
}
