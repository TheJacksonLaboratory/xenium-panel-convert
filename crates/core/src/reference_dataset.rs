use std::fs;

use camino::Utf8Path;
use hdf5_metno::{File, types::FixedAscii};
use ndarray::Array1;
use serde::Serialize;

use crate::{
    error::collect_error,
    reference_dataset::{
        columns::{CellAnnotationCol, CellBarcodeCol, EnsemblIdCol, GeneNameCol},
        error::{
            ReadReferenceDatasetError, ReadReferenceDatasetErrorSet, WriteReferenceDatasetError,
        },
        h5_util::{create_h5_group, write_dataset_to_h5_group},
        obs::{read_cell_annotations_from_h5ad, read_cell_barcodes_from_h5ad},
        pseudo_anndata::PseudoAnndata,
        transcriptome::Transcriptome,
        umi_counts::read_umi_counts_from_h5ad,
        var::read_features_from_h5ad,
    },
};

pub mod columns;
pub mod error;
pub mod h5_util;
pub mod obs;
pub mod pseudo_anndata;
pub mod transcriptome;
pub mod umi_counts;
pub mod var;

pub fn read_reference_dataset(
    path: &Utf8Path,
    cell_barcode_col: &CellBarcodeCol,
    cell_annotation_col: &CellAnnotationCol,
    ensembl_id_col: &EnsemblIdCol,
    gene_name_col: &GeneNameCol,
    transcriptome: Option<Transcriptome>,
) -> Result<PseudoAnndata, ReadReferenceDatasetErrorSet> {
    let mut errors = Vec::new();

    let file = hdf5_metno::File::open(path).map_err(|e| {
        ReadReferenceDatasetErrorSet::new(
            path,
            vec![ReadReferenceDatasetError::H5File {
                reason: e.to_string(),
            }],
        )
    })?;

    let counts = collect_error(read_umi_counts_from_h5ad(&file), &mut errors);

    let barcodes = collect_error(
        read_cell_barcodes_from_h5ad(&file, cell_barcode_col)
            .map_err(ReadReferenceDatasetError::CellBarcodes),
        &mut errors,
    );

    let cell_annotations = collect_error(
        read_cell_annotations_from_h5ad(&file, cell_annotation_col)
            .map_err(ReadReferenceDatasetError::CellAnnotations),
        &mut errors,
    );

    let features = collect_error(
        read_features_from_h5ad(&file, ensembl_id_col, gene_name_col, transcriptome),
        &mut errors,
    );

    let (Some(counts), Some(barcodes), Some(cell_annotations), Some(features)) =
        (counts, barcodes, cell_annotations, features)
    else {
        return Err(ReadReferenceDatasetErrorSet::new(path, errors));
    };

    let anndata = collect_error(
        PseudoAnndata::new(counts, barcodes, cell_annotations, features),
        &mut errors,
    );

    anndata.ok_or_else(|| ReadReferenceDatasetErrorSet::new(path, errors))
}

// See https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices for format specifics
pub fn write_reference_dataset(
    dir: &Utf8Path,
    ds: &PseudoAnndata,
) -> Result<(), WriteReferenceDatasetError> {
    if !dir.exists() {
        fs::create_dir_all(dir).map_err(|e| WriteReferenceDatasetError::CreateOutputDir {
            path: dir.to_owned(),
            reason: e.to_string(),
        })?;
    }

    let annotations_path = dir.join("annotations.csv");
    if annotations_path.exists() {
        return Err(WriteReferenceDatasetError::AnnotationsCsvExists {
            path: annotations_path,
        });
    }
    write_annotations_csv(&annotations_path, ds.barcodes(), ds.cell_annotations())?;

    write_matrix(&dir.join("matrix.h5"), ds)
}

fn write_annotations_csv(
    path: &Utf8Path,
    barcodes: &Barcodes,
    annotations: &CellAnnotations,
) -> Result<(), WriteReferenceDatasetError> {
    #[derive(Debug, Serialize)]
    struct CellAnnotation<'a> {
        barcode: &'a str,
        annotation: &'a str,
    }

    let map_err = |err: csv::Error| WriteReferenceDatasetError::WriteCsv {
        path: path.to_owned(),
        reason: err.to_string(),
    };

    let mut writer = csv::Writer::from_path(path).map_err(map_err)?;
    for (barcode, annotation) in barcodes.iter().zip(annotations) {
        writer
            .serialize(CellAnnotation {
                barcode,
                annotation,
            })
            .map_err(map_err)?;
    }

    Ok(())
}

fn write_matrix(
    path: &Utf8Path,
    dataset: &PseudoAnndata,
) -> Result<(), WriteReferenceDatasetError> {
    let file =
        File::create_excl(path).map_err(|e| WriteReferenceDatasetError::CreateMatrixFile {
            path: path.to_path_buf(),
            reason: e.to_string(),
        })?;

    let write_err = |error| WriteReferenceDatasetError::WriteH5Object {
        path: path.to_path_buf(),
        error,
    };

    let matrix_group = create_h5_group(&file, "matrix").map_err(write_err)?;

    write_dataset_to_h5_group(&matrix_group, "barcodes", dataset.barcodes()).map_err(write_err)?;

    let counts = dataset.counts();
    write_dataset_to_h5_group(&matrix_group, "data", counts.data()).map_err(write_err)?;
    write_dataset_to_h5_group(&matrix_group, "indices", counts.indices()).map_err(write_err)?;
    write_dataset_to_h5_group(&matrix_group, "indptr", counts.indptr().raw_storage())
        .map_err(write_err)?;
    write_dataset_to_h5_group(&matrix_group, "shape", &counts.shape_as_i32()).map_err(write_err)?;

    let features = dataset.features();
    write_dataset_to_h5_group(
        &matrix_group,
        "features/feature_type",
        features.feature_types(),
    )
    .map_err(write_err)?;
    write_dataset_to_h5_group(&matrix_group, "features/id", features.ensembl_ids())
        .map_err(write_err)?;
    write_dataset_to_h5_group(&matrix_group, "features/name", features.gene_names())
        .map_err(write_err)?;

    Ok(())
}

type Barcode = FixedAscii<64>;

type Barcodes = Array1<Barcode>;

type CellAnnotations = Array1<String>;

#[cfg(test)]
mod tests {
    use std::fs;

    use camino::Utf8Path;
    use hdf5_metno::{
        File,
        types::{FixedAscii, TypeDescriptor},
    };
    use tempfile::TempDir;

    use crate::{
        error::Hinted,
        reference_dataset::{
            Barcode,
            columns::{CellAnnotationCol, CellBarcodeCol, EnsemblIdCol, GeneNameCol},
            error::{
                ReadReferenceDatasetError, ReadReferenceDatasetErrorSet, WriteReferenceDatasetError,
            },
            h5_util::read_test_1d_dataset,
            pseudo_anndata::PseudoAnndata,
            read_reference_dataset,
            transcriptome::{Transcriptome, TranscriptomeName},
            var::{EnsemblId, GeneName, VarError},
            write_reference_dataset,
        },
    };

    const REAL_H5AD: &str = "test-data/1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix.h5ad";
    const REAL_H5: &str = "test-data/1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix.h5";

    fn read_scanpy_generated_dataset() -> PseudoAnndata {
        read_reference_dataset(
            Utf8Path::new(REAL_H5AD),
            &CellBarcodeCol("barcode".to_owned()),
            &CellAnnotationCol("annotation".to_owned()),
            &EnsemblIdCol("gene_ids".to_owned()),
            &GeneNameCol("gene_name".to_owned()),
            Transcriptome::new(TranscriptomeName::Mm102020A, false),
        )
        .unwrap()
    }

    fn temp_dir() -> TempDir {
        TempDir::new().unwrap()
    }

    fn utf8_path_from_temp_dir(dir: &TempDir) -> &Utf8Path {
        Utf8Path::from_path(dir.path()).unwrap()
    }

    fn dtype(file: &File, path: &str) -> TypeDescriptor {
        file.dataset(path)
            .unwrap()
            .dtype()
            .unwrap()
            .to_descriptor()
            .unwrap()
    }

    #[test]
    fn read_collects_errors_from_every_field() {
        let dataset_path = "test-data/csr_adata.h5ad";

        let ReadReferenceDatasetErrorSet { path, errors } = read_reference_dataset(
            Utf8Path::new(dataset_path),
            &CellBarcodeCol("foo".to_owned()),
            &CellAnnotationCol("bar".to_owned()),
            &EnsemblIdCol("baz".to_owned()),
            &GeneNameCol("qux".to_owned()),
            Transcriptome::new(TranscriptomeName::Grch382020A, false),
        )
        .unwrap_err();

        assert_eq!(path, dataset_path);

        std::assert_matches!(
            errors.as_slice(),
            [
                Hinted {
                    error: ReadReferenceDatasetError::CellBarcodes(..),
                    ..
                },
                Hinted {
                    error: ReadReferenceDatasetError::CellAnnotations(..),
                    ..
                },
                Hinted {
                    error: ReadReferenceDatasetError::Var(VarError::InvalidH5Fields { .. }),
                    ..
                },
            ],
            "every unreadable field should be reported, not just the first"
        );
    }

    #[test]
    fn written_dataset_is_the_same_as_cellranger_h5() {
        let scanpy_dataset = read_scanpy_generated_dataset();

        let dir = temp_dir();
        let output_dir = utf8_path_from_temp_dir(&dir);
        write_reference_dataset(output_dir, &scanpy_dataset).unwrap();

        let written = File::open(output_dir.join("matrix.h5")).unwrap();

        let scanpy_counts = scanpy_dataset.counts();

        assert_eq!(
            read_test_1d_dataset::<i32>(&written, "matrix/data")
                .unwrap()
                .as_slice()
                .unwrap(),
            scanpy_counts.data()
        );
        assert_eq!(
            read_test_1d_dataset::<i64>(&written, "matrix/indices")
                .unwrap()
                .as_slice()
                .unwrap(),
            scanpy_counts.indices()
        );
        assert_eq!(
            read_test_1d_dataset::<i64>(&written, "matrix/indptr")
                .unwrap()
                .as_slice()
                .unwrap(),
            scanpy_counts.indptr().raw_storage()
        );
        assert_eq!(
            read_test_1d_dataset::<i32>(&written, "matrix/shape")
                .unwrap()
                .as_slice()
                .unwrap(),
            scanpy_counts.shape_as_i32()
        );

        assert_eq!(
            read_test_1d_dataset::<Barcode>(&written, "matrix/barcodes").unwrap(),
            scanpy_dataset.barcodes()
        );

        let scanpy_features = scanpy_dataset.features();
        assert_eq!(
            read_test_1d_dataset::<EnsemblId>(&written, "matrix/features/id").unwrap(),
            scanpy_features.ensembl_ids()
        );
        assert_eq!(
            read_test_1d_dataset::<GeneName>(&written, "matrix/features/name").unwrap(),
            scanpy_features.gene_names()
        );
        assert_eq!(
            read_test_1d_dataset::<FixedAscii<32>>(&written, "matrix/features/feature_type")
                .unwrap(),
            scanpy_features.feature_types()
        );

        // Also check the data types
        let original = File::open(REAL_H5).unwrap();
        for path in [
            "matrix/data",
            "matrix/indices",
            "matrix/indptr",
            "matrix/shape",
        ] {
            assert_eq!(
                dtype(&written, path),
                dtype(&original, path),
                "{path} was not written with the type cellranger writes"
            );
        }

        // Check the annotations CSV too
        let annotations = fs::read_to_string(output_dir.join("annotations.csv")).unwrap();
        let mut rows = annotations.lines();

        assert_eq!(rows.next().unwrap(), "barcode,annotation");

        let rows: Vec<_> = rows.collect();
        assert_eq!(
            rows.len(),
            scanpy_dataset.barcodes().len(),
            "expected one row per barcode"
        );

        // The annotations of the test-data alternate between two groups
        assert!(rows[0].ends_with(",group_0"));
        assert!(rows[1].ends_with(",group_1"));
    }

    #[test]
    fn does_not_overwrite_existing_outputs() {
        let dataset = read_scanpy_generated_dataset();
        let dir = temp_dir();

        let existing_annotations_path = utf8_path_from_temp_dir(&dir).join("existing-annotations");
        fs::create_dir_all(&existing_annotations_path).unwrap();
        fs::write(existing_annotations_path.join("annotations.csv"), "").unwrap();

        std::assert_matches!(
            write_reference_dataset(&existing_annotations_path, &dataset).unwrap_err(),
            WriteReferenceDatasetError::AnnotationsCsvExists { .. }
        );

        let existing_matrix = utf8_path_from_temp_dir(&dir).join("existing-matrix");
        fs::create_dir_all(&existing_matrix).unwrap();
        fs::write(existing_matrix.join("matrix.h5"), "").unwrap();

        std::assert_matches!(
            write_reference_dataset(&existing_matrix, &dataset).unwrap_err(),
            WriteReferenceDatasetError::CreateMatrixFile { .. }
        );
    }

    #[test]
    // "read" in the past-tense, as in "the dataset was read" (pronounced like
    // "red")
    fn read_dataset_is_the_same_as_cellranger_h5() {
        let read_dataset = read_scanpy_generated_dataset();

        // Compare it against the original data
        let original_h5 = hdf5_metno::File::open(REAL_H5).unwrap();

        let read_counts = read_dataset.counts();
        let original_counts = read_test_1d_dataset::<i32>(&original_h5, "matrix/data").unwrap();
        assert_eq!(
            original_counts.as_slice().unwrap(),
            read_counts.data(),
            "UMI counts were not correctly reconstructed"
        );

        let original_indices = read_test_1d_dataset::<i64>(&original_h5, "matrix/indices").unwrap();
        assert_eq!(
            original_indices.as_slice().unwrap(),
            read_counts.indices(),
            "UMI count indices were not correctly? reconstructed"
        );

        let original_indptr = read_test_1d_dataset::<i64>(&original_h5, "matrix/indptr").unwrap();
        assert_eq!(
            original_indptr.as_slice().unwrap(),
            read_counts.indptr().raw_storage(),
            "UMI counts indptr was not correctly? reconstructed"
        );
        assert_eq!(
            read_test_1d_dataset::<i32>(&original_h5, "matrix/shape")
                .unwrap()
                .as_slice()
                .unwrap(),
            read_counts.shape_as_i32()
        );

        let original_barcodes =
            read_test_1d_dataset::<FixedAscii<64>>(&original_h5, "matrix/barcodes").unwrap();
        assert_eq!(original_barcodes, read_dataset.barcodes());

        let read_features = read_dataset.features();
        let original_feature_ids =
            read_test_1d_dataset::<EnsemblId>(&original_h5, "matrix/features/id").unwrap();
        assert_eq!(original_feature_ids, read_features.ensembl_ids());

        let original_feature_names =
            read_test_1d_dataset::<GeneName>(&original_h5, "matrix/features/name").unwrap();
        assert_eq!(original_feature_names, read_features.gene_names());
    }
}
