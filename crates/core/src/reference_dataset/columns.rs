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
pub struct GeneSymbolCol(pub String);

impl GeneSymbolCol {
    #[must_use]
    pub(super) fn as_str(&self) -> &str {
        &self.0
    }
}

impl Display for GeneSymbolCol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl Default for GeneSymbolCol {
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

    #[must_use]
    pub fn new(layer_name: &str) -> Self {
        if layer_name == "X" {
            Self::x()
        } else {
            Self(format!("layers/{layer_name}"))
        }
    }

    #[must_use]
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

#[cfg(test)]
mod tests {
    use crate::reference_dataset::columns::CountsLayerName;

    #[test]
    fn counts_layer_prepends_layers_scope() {
        assert_eq!(CountsLayerName::new("counts").as_str(), "layers/counts");
    }

    #[test]
    fn counts_layer_name_does_not_prepend_layers_scope_if_x() {
        assert_eq!(CountsLayerName::new("X").as_str(), "X");

        assert_eq!(CountsLayerName::x(), CountsLayerName::new("X"));
    }
}
