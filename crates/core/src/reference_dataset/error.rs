use camino::{Utf8Path, Utf8PathBuf};
use serde::Serialize;

use crate::{
    error::Hinted,
    reference_dataset::{
        h5_util::{ReadH5FieldError, WriteH5ObjectError},
        pseudo_anndata::ShapeMismatchError,
        umi_counts::UmiCountsError,
        var::VarError,
    },
};

#[derive(Clone, Debug, Serialize)]
pub struct ReadReferenceDatasetErrorSet {
    pub path: Utf8PathBuf,
    pub errors: Vec<Hinted<ReadReferenceDatasetError>>,
}

impl ReadReferenceDatasetErrorSet {
    pub(super) fn new(path: &Utf8Path, errors: Vec<ReadReferenceDatasetError>) -> Self {
        Self {
            path: path.to_owned(),
            errors: errors.into_iter().map(Hinted::new).collect(),
        }
    }
}

#[derive(Clone, Debug, Serialize, thiserror::Error)]
#[serde(tag = "component", rename_all = "snake_case")]
pub enum ReadReferenceDatasetError {
    #[error(
        "the file could not be opened as an H5 file ({reason}) - ensure it exists and was written \
         by scanpy"
    )]
    H5File { reason: String },
    #[error(transparent)]
    UmiCounts(#[from] UmiCountsError),
    #[error(transparent)]
    CellBarcodes(ReadH5FieldError),
    #[error(transparent)]
    CellAnnotations(ReadH5FieldError),
    #[error(transparent)]
    Var(#[from] VarError),
    #[error(transparent)]
    Shape(#[from] ShapeMismatchError),
}

#[derive(Clone, Debug, Serialize, thiserror::Error)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum WriteReferenceDatasetError {
    #[error("the output directory {path} could not be created ({reason})")]
    CreateOutputDir { path: Utf8PathBuf, reason: String },
    #[error("{path} could not be created ({reason})")]
    CreateMatrixFile { path: Utf8PathBuf, reason: String },
    #[error(
        "the H5 object {} could not be written to file {path} ({})", error.object_path, error.reason
    )]
    WriteH5Object {
        path: Utf8PathBuf,
        error: WriteH5ObjectError,
    },
    #[error("cannot overwrite {path} - move or delete the existing annotations.csv file")]
    AnnotationsCsvExists { path: Utf8PathBuf },
    #[error("cannot write CSV to {path} ({reason})")]
    WriteCsv { path: Utf8PathBuf, reason: String },
}
