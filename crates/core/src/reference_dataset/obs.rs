use hdf5_metno::File;

use crate::reference_dataset::{
    Barcodes, CellAnnotations,
    columns::{CellAnnotationCol, CellBarcodeCol},
    h5_util::{ReadH5FieldError, read_1d_string_dataset, to_ascii},
};

// Read these into String because we will write them to a CSV, so we need serde
// support
pub(super) fn read_cell_annotations_from_h5ad(
    file: &File,
    annotation_col: &CellAnnotationCol,
) -> Result<CellAnnotations, ReadH5FieldError> {
    let strings = read_1d_string_dataset(file, &format!("obs/{annotation_col}"))?;

    Ok(strings.mapv_into_any(|s| s.to_string()))
}

// If the dataset combines multiple smaller datasets, then the barcodes may
// contain sample names, which can be arbitrarily long. As such, we allow a
// conservative 64 bytes per barcode (18 bytes for the barcode itself and 46 for
// the sample name). This is easy to adjust should we find it to be too small
pub(super) fn read_cell_barcodes_from_h5ad(
    file: &File,
    barcode_col: &CellBarcodeCol,
) -> Result<Barcodes, ReadH5FieldError> {
    let barcodes = read_1d_string_dataset(file, &format!("obs/{barcode_col}"))?;

    Ok(barcodes.mapv_into_any(|b| to_ascii(&b)))
}

#[cfg(test)]
mod tests {
    use hdf5_metno::{File, types::FixedAscii};
    use ndarray::Array1;

    use crate::reference_dataset::{
        columns::{CellAnnotationCol, CellBarcodeCol},
        h5_util::{FieldType, ReadH5FieldError},
        obs::{read_cell_annotations_from_h5ad, read_cell_barcodes_from_h5ad},
    };

    fn generated_h5ad() -> File {
        File::open("test-data/csr_adata.h5ad").unwrap()
    }

    #[test]
    fn reads_barcodes_and_annotations() {
        let file = generated_h5ad();

        let barcodes =
            read_cell_barcodes_from_h5ad(&file, &CellBarcodeCol("barcode".to_owned())).unwrap();
        let expected_barcodes: Array1<FixedAscii<64>> = (0..10)
            .map(|i| FixedAscii::from_ascii(&format!("cell_{i}")).unwrap())
            .collect();
        assert_eq!(barcodes, expected_barcodes);

        let annotations =
            read_cell_annotations_from_h5ad(&file, &CellAnnotationCol("annotation".to_owned()))
                .unwrap();
        let expected_annotations: Vec<_> = (0..10).map(|i| format!("group{}", i % 2)).collect();
        assert_eq!(annotations.to_vec(), expected_annotations);
    }

    #[test]
    fn default_barcode_column_reads_the_anndata_index() {
        let file = generated_h5ad();

        let from_index = read_cell_barcodes_from_h5ad(&file, &CellBarcodeCol::default()).unwrap();
        let from_column =
            read_cell_barcodes_from_h5ad(&file, &CellBarcodeCol("barcode".to_owned())).unwrap();

        assert_eq!(
            from_index, from_column,
            "the default barcode column should read anndata's index"
        );
    }

    #[test]
    fn missing_column_is_an_error() {
        let err = read_cell_annotations_from_h5ad(
            &generated_h5ad(),
            &CellAnnotationCol("nonexistent".to_owned()),
        )
        .unwrap_err();

        std::assert_matches!(
            err,
            ReadH5FieldError::InvalidH5ObjectPath {
                field_type: FieldType::Container,
                ..
            }
        );
    }
}
