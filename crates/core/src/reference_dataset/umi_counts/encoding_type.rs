use std::str::FromStr;

#[derive(Clone, Copy, Debug, strum::Display)]
pub(super) enum EncodingType {
    #[strum(transparent)]
    Dense(DenseEncodingType),
    #[strum(transparent)]
    Sparse(SparseEncodingType),
}

impl EncodingType {
    pub(super) const VARIANTS: &'static [&'static str] = &["array", "csr_matrix", "csc_matrix"];
}

// We don't actually need a good error here because the caller decides what to
// do
impl FromStr for EncodingType {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match SparseEncodingType::from_str(s).map(Self::Sparse) {
            Ok(enc) => Ok(enc),
            Err(_) => DenseEncodingType::from_str(s).map(Self::Dense).map_err(|_| ()),
        }
    }
}

#[derive(Clone, Copy, Debug, strum::EnumString, strum::Display)]
#[strum(serialize_all = "snake_case")]
pub(super) enum DenseEncodingType {
    Array,
}

#[derive(Clone, Copy, Debug, strum::EnumString, strum::Display)]
#[strum(serialize_all = "snake_case")]
pub(super) enum SparseEncodingType {
    CsrMatrix,
    CscMatrix,
}
