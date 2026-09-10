use std::{
    collections::{HashMap, HashSet},
    str::FromStr,
};

use anyhow::{anyhow, bail, ensure};
use camino::{Utf8Path, Utf8PathBuf};
use xenium_panel_convert_core::reference_dataset::{
    columns::{CellAnnotationCol, CellBarcodeCol, CountsLayerName, EnsemblIdCol, GeneSymbolCol},
    pseudo_anndata::PseudoAnndata,
    read_reference_dataset,
    transcriptome::{Transcriptome, TranscriptomeName},
    write_reference_dataset,
};

use crate::write::write_json_to_file;

pub(super) fn convert_all_reference_datasets_and_write<'a>(
    ReferenceDatasetCliOptions { reference_datasets }: &'a ReferenceDatasetCliOptions,
    output_dir: &Utf8Path,
) -> anyhow::Result<()> {
    for spec in reference_datasets {
        convert_reference_dataset_and_write(spec, output_dir)?;
    }

    Ok(())
}

pub(super) fn convert_reference_dataset_and_write<'a>(
    spec: &'a ReferenceDatasetSpec,
    output_dir: &Utf8Path,
) -> anyhow::Result<PseudoAnndata> {
    let ReferenceDatasetSpec {
        path,
        counts_layer_name,
        cell_barcode_col,
        cell_annotation_col,
        ensembl_id_col,
        gene_symbol_col,
        transcriptome_name: _,
        flex: _,
        transcriptome,
        rename,
    } = spec;
    let dataset_name = dataset_name(path, rename.as_deref())?;

    let output_path = |filename: &str| output_dir.join(filename);

    match read_reference_dataset(
        path,
        counts_layer_name,
        cell_barcode_col,
        cell_annotation_col,
        ensembl_id_col,
        gene_symbol_col,
        *transcriptome,
    ) {
        Ok(ds) => {
            write_reference_dataset(&output_path(dataset_name), &ds)?;
            Ok(ds)
        }
        Err(errors) => {
            let error_path = format!("{dataset_name}-errors.json");
            write_json_to_file(&errors, &output_path(&error_path))?;

            Err(errors.into())
        }
    }
}

pub(super) fn dataset_name<'a>(path: &'a Utf8Path, rename: Option<&'a Utf8Path>) -> anyhow::Result<&'a str> {
    rename
        .and_then(Utf8Path::file_stem)
        .or_else(|| path.file_stem())
        .ok_or_else(|| anyhow!("invalid filename: {}", rename.unwrap_or(path)))
}

#[derive(Clone, Debug)]
pub(super) struct ReferenceDatasetSpec {
    pub(super) path: Utf8PathBuf,
    counts_layer_name: CountsLayerName,
    cell_barcode_col: CellBarcodeCol,
    cell_annotation_col: CellAnnotationCol,
    ensembl_id_col: EnsemblIdCol,
    gene_symbol_col: GeneSymbolCol,
    pub(super) transcriptome_name: TranscriptomeName,
    pub(super) flex: bool,
    transcriptome: Option<Transcriptome>,
    pub(super) rename: Option<Utf8PathBuf>,
}

impl ReferenceDatasetSpec {
    fn parse_commandline(s: &str) -> anyhow::Result<Self> {
        const EXAMPLE: &str = "path=matrix.h5ad,counts-layer=X,barcode-col=barcodes,annotation-col=annotations,\
                               ensembl-id-col=gene_ids,transcriptome=GRCh38-2024-A\nmatrix.h5ad,b=barcodes,\
                               a=annotations,e=gene_ids,t=h2024";

        fn get_spec_value_default<T: Default>(spec: &HashMap<&str, &str>, key: &str, map: impl Fn(String) -> T) -> T {
            spec.get(key).map(|&s| s.to_owned()).map(map).unwrap_or_default()
        }

        fn get_spec_value<T>(spec: &HashMap<&str, &str>, key: &str, map: impl Fn(String) -> T) -> anyhow::Result<T> {
            let val = *spec.get(key).ok_or_else(|| anyhow!("key '{key}' is required"))?;

            Ok(map(val.to_owned()))
        }

        let key_aliases: HashMap<_, _> = [
            ("p", "path"),
            ("c", "counts-layer"),
            ("b", "barcode-col"),
            ("a", "annotation-col"),
            ("e", "ensembl-id-col"),
            ("g", "gene-symbol-col"),
            ("t", "transcriptome"),
            ("f", "flex"),
            ("r", "rename"),
        ]
        .into_iter()
        .collect();
        let allowed_keys: HashSet<_> = key_aliases.iter().flat_map(|(s1, s2)| [*s1, *s2]).collect();

        let mut spec = HashMap::with_capacity(6);

        for (i, kv_pair) in s.split(',').enumerate() {
            let (key, value) = match (i, kv_pair.split_once('=')) {
                (_, Some((key, value))) => (key, value),
                (0, None) => ("path", kv_pair),
                (_, None) => {
                    bail!(
                        "reference dataset specification must be provided like one of the following (mixing aliases \
                         and full-names is allowed):\n{EXAMPLE}"
                    );
                }
            };

            ensure!(
                allowed_keys.contains(key),
                "key '{key}' not recognized in reference dataset specification"
            );

            let key = key_aliases.get(key).unwrap_or(&key);

            ensure!(
                spec.insert(*key, value).is_none(),
                "key '{key}' may not be specified more than once"
            );
        }

        let transcriptome_from_str = |s: &&str| TranscriptomeName::from_str(s).map_err(anyhow::Error::from);
        let transcriptome_name = spec
            .get("transcriptome")
            .ok_or(anyhow!("key 'transcriptome' is required"))
            .and_then(transcriptome_from_str)?;
        let flex = spec
            .get("flex")
            .map(|s| bool::from_str(s))
            .transpose()?
            .unwrap_or_default();

        Ok(Self {
            path: get_spec_value(&spec, "path", Utf8PathBuf::from)?,
            counts_layer_name: get_spec_value_default(&spec, "counts-layer", |s| CountsLayerName::new(&s)),
            cell_barcode_col: get_spec_value_default(&spec, "barcode-col", CellBarcodeCol),
            cell_annotation_col: get_spec_value(&spec, "annotation-col", CellAnnotationCol)?,
            ensembl_id_col: get_spec_value_default(&spec, "ensembl-id-col", EnsemblIdCol),
            gene_symbol_col: get_spec_value_default(&spec, "gene-symbol-col", GeneSymbolCol),
            transcriptome_name,
            flex,
            transcriptome: Transcriptome::new(transcriptome_name, flex),
            rename: get_spec_value_default(&spec, "rename", |s| Some(Utf8PathBuf::from(s))),
        })
    }
}

#[derive(Clone, Debug, clap::Args)]
pub(crate) struct ReferenceDatasetCliOptions {
    #[clap(
        value_parser = ReferenceDatasetSpec::parse_commandline,
        help = REFERENCE_DATASETS_HELP,
        long_help = REFERENCE_DATASETS_LONG_HELP
    )]
    pub(super) reference_datasets: Vec<ReferenceDatasetSpec>,
}

const REFERENCE_DATASETS_HELP: &str =
    "One or more reference datasets, each specified as a comma-delimited list of <KEY>=<VALUE> pairs.";

const REFERENCE_DATASETS_LONG_HELP: &str =
    "One or more reference datasets, each specified as a comma-delimited list of <KEY>=<VALUE> pairs. Every key has a \
     one-letter alias, and full names and aliases may be mixed within a single specification. No key may be given \
     more than once. The value of 'path' may also be given without its key, in which case it must come first.

KEY              ALIAS  VALUE
- - - - - - - - - - - - - - -
path            | p | the h5ad file to convert (required)
annotation-col  | a | obs column containing cell annotations (required)
transcriptome   | t | transcriptome against which the dataset was aligned (required)
counts-layer    | c | layer where counts are stored [default: X]
barcode-col     | b | obs column containing cell barcodes [default: _index]
ensembl-id-col  | e | var column containing Ensembl IDs [default: gene_ids]
gene-symbol-col | g | var column containing gene symbols [default: _index]
flex            | f | 'true' if the dataset was generated with the Flex (probe-based) chemistry, 'false' otherwise \
     [default: false]
rename          | r | name of the converted dataset in <OUTPUT_DIR> [default: the filename of 'path', without its \
     extension]


'transcriptome' is case-insensitive and must be one of:
NAME           ALIAS
- - - - - - - - - - -
GRCh38-2020-A | h2020
GRCh38-2024-A | h2024
mm10-2020-A   | m2020
GRCm39-2024-A | m2024
other         | o

Note that setting 'transcriptome' = 'other' will disable checking .var for gene-filtering.

Examples:

- path=matrix.h5ad,annotation-col=cell_type,transcriptome=GRCh38-2024-A
- matrix.h5ad,a=cell_type,t=h2024,f=true,r=renamed-dataset";

#[cfg(test)]
mod tests {
    use camino::Utf8Path;
    use xenium_panel_convert_core::reference_dataset::{
        columns::{CellAnnotationCol, CellBarcodeCol, EnsemblIdCol, GeneSymbolCol},
        transcriptome::Transcriptome,
    };

    use crate::reference_datasets::{ReferenceDatasetSpec, dataset_name};

    #[test]
    fn spec_parses_positional_path_and_aliases() {
        let spec = ReferenceDatasetSpec::parse_commandline(
            "matrix.h5ad,b=barcode,a=annotation,e=id,g=name,t=h2020,f=true,r=renamed",
        )
        .unwrap();

        assert_eq!(spec.path.as_str(), "matrix.h5ad");
        assert_eq!(spec.cell_barcode_col, CellBarcodeCol("barcode".to_owned()));
        assert_eq!(spec.cell_annotation_col, CellAnnotationCol("annotation".to_owned()));
        assert_eq!(spec.ensembl_id_col, EnsemblIdCol("id".to_owned()));
        assert_eq!(spec.gene_symbol_col, GeneSymbolCol("name".to_owned()));
        assert_eq!(spec.rename.as_deref(), Some(Utf8Path::new("renamed")));
        std::assert_matches!(spec.transcriptome.unwrap(), Transcriptome::Flex2020A(_));
    }

    #[test]
    fn spec_applies_column_defaults() {
        let spec =
            ReferenceDatasetSpec::parse_commandline("path=matrix.h5ad,annotation-col=annotations,transcriptome=m2020")
                .unwrap();

        assert_eq!(spec.cell_barcode_col, CellBarcodeCol("_index".to_owned()));
        assert_eq!(spec.ensembl_id_col, EnsemblIdCol("gene_ids".to_owned()));
        assert_eq!(spec.gene_symbol_col, GeneSymbolCol("_index".to_owned()));
        assert_eq!(spec.rename, None);
    }

    #[test]
    fn spec_rejects_malformed_specs() {
        let malformed = [
            "matrix.h5ad,a=annotations,t=h2020,unknown=foo",
            "p=matrix.h5ad,path=other.h5ad,a=annotations,t=h2020",
            "matrix.h5ad,a=annotations,t=h2020,bare-value",
            "matrix.h5ad,t=h2020",
            "matrix.h5ad,a=annotations,t=not-a-transcriptome",
        ];

        for spec in malformed {
            assert!(
                ReferenceDatasetSpec::parse_commandline(spec).is_err(),
                "'{spec}' should not have parsed"
            );
        }
    }

    #[test]
    fn dataset_names_come_from_file_stems() {
        let path = Utf8Path::new("datasets/matrix.h5ad");

        assert_eq!(dataset_name(path, None).unwrap(), "matrix");
        assert_eq!(
            dataset_name(path, Some(Utf8Path::new("elsewhere/renamed.h5ad"))).unwrap(),
            "renamed",
            "rename should take precedence over the dataset's path"
        );
    }
}
