use camino::Utf8PathBuf;
use serde::Serialize;

use crate::{
    common::ErrorVecExt,
    reference_dataset::{
        h5_util::{CreateH5GroupError, ReadH5FieldError, WriteH5DatasetError},
        pseudo_anndata::ShapeMismatchError,
        umi_counts::UmiCountsError,
        var::VarError,
    },
};

#[derive(Clone, Debug, Serialize, thiserror::Error)]
#[serde(untagged)]
pub enum ReadReferenceDatasetErrorSet {
    #[error("invalid H5 file at {path} ({reason}) - {hint}")]
    InvalidH5File {
        path: Utf8PathBuf,
        reason: String,
        hint: &'static str,
    },
    #[error("invalid matrix")]
    Matrix {
        path: Utf8PathBuf,
        errors: Vec<ReadReferenceDatasetError>,
    },
}

#[derive(Clone, Debug, Serialize, thiserror::Error)]
#[error("{error}")]
pub struct ReadReferenceDatasetError {
    #[serde(flatten)]
    pub error: ReadReferenceDatasetErrorInner,
    #[serde(skip_serializing_if = "str::is_empty")]
    pub hint: String,
}

#[derive(Clone, Debug, Serialize, thiserror::Error)]
#[serde(tag = "component", rename_all = "snake_case")]
pub enum ReadReferenceDatasetErrorInner {
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

impl<E> ErrorVecExt<E> for Vec<ReadReferenceDatasetErrorInner>
where
    E: Into<ReadReferenceDatasetErrorInner>,
{
    fn push_err<T>(&mut self, err: E) -> Option<T> {
        self.push(err.into());

        None
    }
}

#[derive(Clone, Debug, Serialize, thiserror::Error)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum WriteReferenceDatasetErrorInner {
    #[error("failed to create output directory - {reason}")]
    CreateOutputDir { path: Utf8PathBuf, reason: String },
    #[error("failed to create matrix.h5 - {reason}")]
    CreateMatrixFile { path: Utf8PathBuf, reason: String },
    #[error("{error}")]
    CreateH5Group {
        path: Utf8PathBuf,
        error: CreateH5GroupError,
    },
    #[error("{error}")]
    WriteH5Dataset {
        path: Utf8PathBuf,
        error: WriteH5DatasetError,
    },
    #[error("cannot overwrite {path} - move or delete the existing annotation.csv file")]
    AnnotationsCsvExists { path: Utf8PathBuf },
}

impl From<ReadReferenceDatasetErrorInner> for ReadReferenceDatasetError {
    fn from(error: ReadReferenceDatasetErrorInner) -> Self {
        Self {
            hint: error.to_string(),
            error,
        }
    }
}

#[derive(Clone, Debug, Serialize, thiserror::Error)]
#[error("{error}")]
pub struct WriteReferenceDatasetError {
    pub error: WriteReferenceDatasetErrorInner,
    pub hint: String,
}

impl From<WriteReferenceDatasetErrorInner> for WriteReferenceDatasetError {
    fn from(error: WriteReferenceDatasetErrorInner) -> Self {
        Self {
            hint: error.to_string(),
            error,
        }
    }
}
