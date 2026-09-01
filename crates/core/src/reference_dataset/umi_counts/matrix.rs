use std::ops::Range;

use ndarray::Array2;
use sprs::{CsMatI, IndPtrView};

use crate::reference_dataset::umi_counts::{
    UmiCountsError,
    encoding_type::{DenseEncodingType, SparseEncodingType},
    matrix::csc::CscMatrix,
};

mod csc;

type Matrix<N> = CsMatI<N, i64>;

type UnvalidatedCscMatrix = CscMatrix<f32>;

type ValidatedCscMatrix = CscMatrix<i32>;

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct RawCscUmiCounts(ValidatedCscMatrix);

impl RawCscUmiCounts {
    pub(super) fn from_sparse_matrix(
        shape: (usize, usize),
        indptr: Vec<i64>,
        indices: Vec<i64>,
        data: Vec<f32>,
        encoding_type: SparseEncodingType,
    ) -> Result<Self, UmiCountsError> {
        let mtx = match encoding_type {
            SparseEncodingType::CsrMatrix => Matrix::try_new(shape, indptr, indices, data),
            SparseEncodingType::CscMatrix => Matrix::try_new_csc(shape, indptr, indices, data),
        };

        mtx.map(Matrix::transpose_into)
            .map(CscMatrix::new)
            .map_err(|(.., err)| err.into())
            .and_then(Self::from_csc_matrix)
    }

    fn from_csc_matrix(mtx: UnvalidatedCscMatrix) -> Result<Self, UmiCountsError> {
        let mtx = validate_counts(mtx)?;
        Ok(Self(mtx))
    }

    pub(super) fn from_dense_matrix(
        counts: &Array2<f32>,
        encoding_type: DenseEncodingType,
    ) -> Result<Self, UmiCountsError> {
        match encoding_type {
            DenseEncodingType::Array => Self::from_csc_matrix(CscMatrix::new(
                Matrix::csc_from_dense(counts.view(), 0.).transpose_into(),
            )),
        }
    }

    pub(crate) fn data(&self) -> &[i32] {
        self.0.as_matrix().data()
    }

    pub(crate) fn indices(&self) -> &[i64] {
        self.0.as_matrix().indices()
    }

    pub(crate) fn indptr(&self) -> IndPtrView<'_, i64> {
        self.0.as_matrix().indptr()
    }

    pub(crate) fn shape(&self) -> [i32; 2] {
        let (nrows, ncols) = self.0.as_matrix().shape();

        #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
        [nrows as i32, ncols as i32]
    }
}

fn validate_counts(mtx: UnvalidatedCscMatrix) -> Result<ValidatedCscMatrix, UmiCountsError> {
    let inner_mtx = mtx.as_matrix();
    let shape = inner_mtx.shape();

    let i32_data: Vec<_> = inner_mtx
        .data()
        .iter()
        .map(|f| f32_to_i32(*f))
        .collect::<Result<_, _>>()?;

    if all_total_counts_are_equal(&i32_data, inner_mtx.indptr())? {
        return Err(UmiCountsError::NormalizedCounts);
    }

    let (indptr, indices, _) = mtx.into_inner().into_raw_storage();
    Ok(CscMatrix::new(Matrix::new_csc(
        shape, indptr, indices, i32_data,
    )))
}

fn all_total_counts_are_equal(
    data: &[i32],
    indptr: IndPtrView<'_, i64>,
) -> Result<bool, UmiCountsError> {
    let total_counts = calculate_total_counts_for_all_cells(data, indptr);

    let mut nonempty_cells = total_counts.iter().filter(|tc| **tc > 0);
    let first_nonempty_cell = nonempty_cells.next().ok_or(UmiCountsError::EmptyCounts)?;

    Ok(nonempty_cells.all(|tc| tc == first_nonempty_cell))
}

fn calculate_total_counts_for_all_cells(data: &[i32], indptr: IndPtrView<'_, i64>) -> Vec<i32> {
    let n_cells = indptr.len() - 1;
    let mut total_counts = vec![0; n_cells];

    for (cell_idx, total_counts_for_this_cell) in total_counts.iter_mut().enumerate() {
        let index_range = indptr.outer_inds(cell_idx);
        *total_counts_for_this_cell = calculate_total_counts_for_cell(data, index_range);
    }

    total_counts
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn calculate_total_counts_for_cell(
    data: &[i32],
    Range {
        start: first_count_idx,
        end: last_count_idx,
    }: Range<i64>,
) -> i32 {
    let cell_counts = &data[first_count_idx as usize..last_count_idx as usize];
    cell_counts.iter().sum()
}

#[allow(clippy::float_cmp)]
fn f32_to_i32(f: f32) -> Result<i32, UmiCountsError> {
    let is_nonnegative_integral = f.round() == f && f >= 0.0;

    if is_nonnegative_integral {
        #[allow(clippy::cast_possible_truncation)]
        Ok(f as i32)
    } else {
        Err(UmiCountsError::TransformedCounts)
    }
}

#[cfg(test)]
mod tests {
    use ndarray::{Array2, array};

    use crate::reference_dataset::umi_counts::{
        RawCscUmiCounts, UmiCountsError,
        encoding_type::{DenseEncodingType, SparseEncodingType},
        matrix::{Matrix, calculate_total_counts_for_all_cells},
    };

    fn counts() -> Array2<f32> {
        array![[0., 1., 2.], [2., 4., 6.]]
    }

    fn csr() -> Matrix<f32> {
        Matrix::csr_from_dense(counts().view(), 0.)
    }

    fn csc() -> Matrix<f32> {
        Matrix::csc_from_dense(counts().view(), 0.)
    }

    #[test]
    fn total_counts_are_correct() {
        let mtx = csr();
        #[expect(clippy::cast_possible_truncation)]
        let data: Vec<_> = mtx.data().iter().map(|f| *f as i32).collect();

        assert_eq!(
            // This function assumes column-orientation - that is, each row is a feature and each
            // column is a cell
            calculate_total_counts_for_all_cells(&data, mtx.indptr()),
            [3, 12]
        );
    }

    #[allow(clippy::similar_names)]
    #[test]
    fn storage_orders_are_equivalent() {
        let csr = csr();
        let shape = csr.shape();
        let (indptr, indices, data) = csr.into_raw_storage();
        let from_csr = RawCscUmiCounts::from_sparse_matrix(
            shape,
            indptr,
            indices,
            data,
            SparseEncodingType::CsrMatrix,
        )
        .unwrap();

        let csc = csc();
        let shape = csc.shape();
        let (indptr, indices, data) = csc.into_raw_storage();
        let from_csc = RawCscUmiCounts::from_sparse_matrix(
            shape,
            indptr,
            indices,
            data,
            SparseEncodingType::CscMatrix,
        )
        .unwrap();

        assert_eq!(
            from_csr, from_csc,
            "the same counts in CSR and CSC did not produce the same matrix"
        );

        let from_dense =
            RawCscUmiCounts::from_dense_matrix(&counts(), DenseEncodingType::Array).unwrap();
        assert_eq!(
            from_csr, from_dense,
            "dense matrix and sparse matrix provided from the same counts did not produce the \
             same matrix"
        );
    }

    #[test]
    fn equal_cell_sums_return_normalized_error() {
        // Both cells sum to 3
        let normalized_counts = array![[0., 1., 2.], [1., 1., 1.]];

        std::assert_matches!(
            RawCscUmiCounts::from_dense_matrix(&normalized_counts, DenseEncodingType::Array),
            Err(UmiCountsError::NormalizedCounts)
        );
    }

    #[test]
    fn all_zero_counts_returns_empty_counts_error() {
        let all_zero_counts = Array2::zeros((2, 2));

        std::assert_matches!(
            RawCscUmiCounts::from_dense_matrix(&all_zero_counts, DenseEncodingType::Array),
            Err(UmiCountsError::EmptyCounts)
        );
    }

    #[test]
    fn transformed_counts_error() {
        let counts = array![[1., 2.5, 3.], [4., 5., 6.]];

        std::assert_matches!(
            RawCscUmiCounts::from_dense_matrix(&counts, DenseEncodingType::Array),
            Err(UmiCountsError::TransformedCounts)
        );
    }

    #[test]
    fn negative_counts_error() {
        let counts = array![[-1., 2., 3.], [4., 5., 6.]];

        std::assert_matches!(
            RawCscUmiCounts::from_dense_matrix(&counts, DenseEncodingType::Array),
            Err(UmiCountsError::TransformedCounts)
        );
    }
}
