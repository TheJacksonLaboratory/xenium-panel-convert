use std::str::FromStr;

use hdf5_metno::{File, types::VarLenUnicode};
pub(super) use matrix::RawCscUmiCounts;
use serde::Serialize;

use crate::reference_dataset::{
    h5_util::{FieldType, ReadH5FieldError, read_attribute, read_container, read_dataset_raw},
    umi_counts::encoding_type::{DenseEncodingType, EncodingType, SparseEncodingType},
};

mod encoding_type;
mod matrix;

pub(super) fn read_umi_counts_from_h5ad(file: &File) -> Result<RawCscUmiCounts, UmiCountsError> {
    let encoding_type: VarLenUnicode =
        read_attribute(&read_container(file, "X")?, "encoding-type")?;

    let encoding_type = EncodingType::from_str(&encoding_type).map_err(|()| {
        UmiCountsError::UnknownEncodingType {
            found: encoding_type.to_string(),
            expected: EncodingType::VARIANTS,
        }
    })?;

    match encoding_type {
        EncodingType::Sparse(enc) => read_x_sparse(file, enc),
        EncodingType::Dense(enc) => read_x_dense(file, enc),
    }
}

fn read_x_sparse(
    file: &File,
    encoding_type: SparseEncodingType,
) -> Result<RawCscUmiCounts, UmiCountsError> {
    let data = read_dataset_raw(file, "X/data")?;
    let indptr = read_dataset_raw(file, "X/indptr")?;
    let indices = read_dataset_raw(file, "X/indices")?;

    // It's very nice that scanpy decides to store the shape as an attribute rather
    // than following 10x Genomics and storing it as a dataset. It's great when a library built to analyze data changes the format of the data for no discernible reason :)
    let shape = file
        .group("X")
        .and_then(|x| x.attr("shape"))
        .and_then(|sh| sh.read_1d())
        .map_err(|err| {
            ReadH5FieldError::new_invalid_field(err, file, "X/shape", FieldType::Attribute)
        })?;
    let shape = (shape[0], shape[1]);

    RawCscUmiCounts::from_sparse_matrix(shape, indptr, indices, data, encoding_type)
}

fn read_x_dense(
    file: &File,
    encoding_type: DenseEncodingType,
) -> Result<RawCscUmiCounts, UmiCountsError> {
    let counts = file.dataset("X").and_then(|ds| ds.read_2d()).map_err(|_| {
        UmiCountsError::MalformedMatrix {
            reason: "failed to read counts in dataset 'X' as 2D array".to_owned(),
        }
    })?;

    RawCscUmiCounts::from_dense_matrix(&counts, encoding_type)
}

#[derive(Debug, Clone, Serialize, thiserror::Error)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum UmiCountsError {
    #[error(transparent)]
    MalformedCounts { error: ReadH5FieldError },
    #[error(
        "the counts in X have an unknown encoding type {found}, expected one of {expected:?} - \
         ensure the file was written by scanpy"
    )]
    UnknownEncodingType {
        found: String,
        expected: &'static [&'static str],
    },
    #[error("every cell has a total count of zero - provide a dataset with raw UMI counts")]
    EmptyCounts,
    #[error("the counts are not whole, non-negative numbers - provide raw, untransformed counts")]
    TransformedCounts,
    #[error(
        "every cell has the same total count, so the counts have been normalized - provide raw, \
         unnormalized counts"
    )]
    NormalizedCounts,
    #[error("the counts matrix is malformed ({reason}) - ensure the file was written by scanpy")]
    MalformedMatrix { reason: String },
}

impl From<ReadH5FieldError> for UmiCountsError {
    fn from(error: ReadH5FieldError) -> Self {
        Self::MalformedCounts { error }
    }
}

impl From<sprs::errors::StructureError> for UmiCountsError {
    fn from(err: sprs::errors::StructureError) -> Self {
        Self::MalformedMatrix {
            reason: err.to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use hdf5_metno::File;

    use crate::reference_dataset::umi_counts::read_umi_counts_from_h5ad;

    #[test]
    fn read_h5ad_files() {
        // The last one apparently has some compression applied to it because of
        // scanpy.write's default behavior. It's really nice that scanpy's default
        // behavior differs from anndata's default behavior :)
        let files = [
            "csr_adata",
            "csc_adata",
            "dense_adata",
            "1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix",
        ]
        .map(|fname| format!("test-data/{fname}.h5ad"))
        .map(|path| File::open(path).unwrap());

        let mut all_counts = Vec::with_capacity(files.len());

        for f in files {
            let filename = f.filename();
            let counts = read_umi_counts_from_h5ad(&f).unwrap();

            if filename.contains("adata") {
                assert_eq!(
                    counts.data()[0],
                    10,
                    "first entry in UMI counts of {filename} != 10"
                );
            }

            all_counts.push(counts);
        }

        // We know the first 3 files are generated from the same data
        assert_eq!(all_counts[0], all_counts[1]);
        assert_eq!(all_counts[0], all_counts[2]);
    }
}
