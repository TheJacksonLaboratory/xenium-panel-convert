# xenium-panel-convert

A command-line utility for converting files to the formats accepted by the 10x Genomics Xenium Panel Designer.

In order to create a custom Xenium Panel, the [10x Genomics Xenium Panel Designer](https://www.10xgenomics.com/support/software/xenium-panel-designer/latest) requires a list of targets (typically genes) and one or more [reference datasets](https://www.10xgenomics.com/support/software/xenium-panel-designer/latest/tutorials/create-single-cell-reference). This command-line tool takes as input the list of targets as a CSV-file and one or more [scanpy](https://github.com/scverse/scanpy)-generated H5AD files and converts them to the formats accepted by the Xenium Panel Designer.

## Installation

Install the latest version from the [releases page](https://github.com/TheJacksonLaboratory/xenium-panel-convert/releases). The installer script will also install `xenium-panel-convert-update` which updates `xp-convert` to the latest version.

## Usage

`xp-convert` has 3 subcommands:

- [`xp-convert targets`](#xp-convert-targets)
- [`xp-convert references`](#xp-convert-references)
- [`xp-convert all`](#xp-convert-all) **(recommended)**

### **`xp-convert targets`**

#### **Inputs**

The main input for this command is a CSV-formatted target-list. The target-list must contain the following fields:
| Field | Description | Allowed Values | Example | Required |
| --- | --- | --- | --- | --- |
| `ensembl_id` | The Ensembl ID of the target. If the row has `custom == true`, this will not be validated against the reference genome. | Any Ensembl ID from the [allowed list of genes](https://www.10xgenomics.com/support/software/xenium-panel-designer/latest/tutorials/create-gene-list#yesprobe), unless `custom == true`, in which case any string | ENSG00000141510 | yes |
| `gene_symbol` | The gene-symbol of the target. If the row has `custom == true`, this will not be validated against the reference genome. | Any gene-symbol from the [allowed list of genes](https://www.10xgenomics.com/support/software/xenium-panel-designer/latest/tutorials/create-gene-list#yesprobe), unless `custom == true`, in which case any string | TP53 | yes |
| `group` | A key to create various gene "groups". If the Xenium Panel Designer cannot include a target, you can replace it with a target with the same `group`. This can string can be whatever you want, though it may be useful to use biologically-meaningful terms. Note that all group-names will be lowercased in the ouptut. | Any string | senescense | yes |
| `priority` | A priority to assign to the target, used to sort the output and prevent the Xenium Panel Designer from dropping must-have genes from the panel. | `backup` \| `desired` \| `must_have` | `must_have` | yes |
| `custom` | Whether this is a "custom" target, in which case it will not be validated against the reference genome. | `true` \| `false` | `true` | no, defaults to `false` |

Any additional fields will be propagated to the output untouched. This is useful for notes or any information you want to attach to targets.

##### **Example target-list**

```csv
ensembl_id,gene_symbol,group,priority,notes
ENSG00000141510,TP53,tumor,must_have,some interesting fact
```

If the target-list has different fieldnames and you don't want to edit the target-list, you can provide a set of field-aliases that map the fieldname in the file to one of the canonical fieldnames above. For example, if the field containing Ensembl IDs is called "gene ID", you can do:

```bash
xp-convert targets --targets-path /path/to/targets.csv --field-alias 'gene ID=ensembl_id'
```

If you have a common set of aliases, you can store them in a TOML file mapping the alias to the canonical fieldname:

```toml
"gene ID" = "ensembl_id"
```

and then pass that file on the command-line:

```bash
xp-convert targets --targets-path /path/to/targets.csv --field-alias-file /path/to/aliases.toml
```

If the two methods are combined, aliases passed on the command-line take precedence over aliases in the file.

#### **Outputs**

Outputs will be written to the directory provided as the `--output-dir` option, which defaults to `xp-convert`. If there are no errors, it will look like:

```
xp-convert
├── validated-targets.csv
└── xenium-panel-designer-targets.csv
```

`validated-targets.csv` is a copy of the input file with renamed fields and lowercased group-names, while `xenium-panel-designer-targets.csv` is a sorted list of targets suitable for copy-pasting into the panel designer. For example, the input file shown above would result in:

```csv
Gene,Ensembl ID,Probe sets,Force
TP53,ENSG00000141510,,forced
```

All encountered errors are collected and written `<OUTPUT_DIR>/target-list-errors.json`. For example, a row without the field `priority` and with an incorrect gene-symbol would result in:

```json
[
  {
    "line_number": 2,
    "submitted_target": {
      "ensembl_id": "ENSG00000141510",
      "gene_symbol": "SOME GENE SYMBOL",
      "group": "group0",
      "priority": null,
      "custom": null
    },
    "errors": [
      {
        "type": "missing_field",
        "fieldname": "priority",
        "hint": "the field 'priority' is missing - add it to the CSV"
      },
      {
        "type": "ensembl_id_gene_name_mismatch",
        "ensembl_id": "ENSG00000141510",
        "correct_gene_symbol": "TP53",
        "hint": "the gene name corresponding to the Ensembl ID ENSG00000141510 is TP53 - change either the Ensembl ID or the gene name so they match"
      }
    ]
  }
]
```

### **`xp-convert references`**

#### **Inputs**

The main input for this command is one or more H5AD files generated using scanpy. The `AnnData` object must contain cell-barcodes, cell-annotations, Ensembl IDs, and gene-symbols. All of these are automatically populated for you by scanpy besides cell-annotations. You can generate such a file from the outputs of [cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) like so:

```python
import scanpy as sc

adata = sc.read_10x_h5("/path/to/data.h5")

# Whatever logic you want to annotate each cell
adata.obs["annotation"] = ...
```

#### **Outputs**

If no errors are found, the following invocation:

```bash
xp-convert references 1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix.h5ad,annotation-col=annotation,transcriptome=mm10-2020-A
```

would produce the following output:

```
xp-convert
└── 1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix
    ├── annotations.csv
    └── matrix.h5
```

This directory is ready for upload to the Xenium Panel Designer in one of the [accepted formats](https://www.10xgenomics.com/support/software/xenium-panel-designer/latest/tutorials/create-single-cell-reference#ref-formats).

All encountered errors are collected and written to `<OUTPUT_DIR>/<DATASET_NAME>-errors.json`. For example, supplying the wrong column for cell-annotations would result in:

```json
{
  "path": "crates/core/test-data/1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix.h5ad",
  "errors": [
    {
      "component": "cell_annotations",
      "type": "invalid_h5_object_path",
      "hdf5_error": "H5Dopen2(): unable to synchronously open dataset: object 'foo' doesn't exist",
      "object_path": "obs/foo",
      "field_type": "container",
      "available_objects": [
        "/X/data",
        "/X/indices",
        "/X/indptr",
        "/obs/_index/mask",
        "/obs/_index/values",
        "/obs/annotation/categories",
        "/obs/annotation/codes",
        "/obs/barcode/mask",
        "/obs/barcode/values",
        "/var/_index/mask",
        "/var/_index/values",
        "/var/feature_types/categories",
        "/var/feature_types/codes",
        "/var/gene_ids/mask",
        "/var/gene_ids/values",
        "/var/gene_symbol/categories",
        "/var/gene_symbol/codes",
        "/var/genome/categories",
        "/var/genome/codes"
      ],
      "hint": "obs/foo could not be read as a container (H5Dopen2(): unable to synchronously open dataset: object 'foo' doesn't exist) - ensure the correct column name was provided"
    }
  ]
}
```

### **`xp-convert all`**

This command is the same as combining the above commands with the added checks that:

- the species of the target-list and reference datasets match
- the gene-symbols in the target-list and reference datasets match
- all genes in the target-list are in the reference datasets

If any of these conditions are not met, warnings are written to `<OUTPUT_DIR>/<DATASET_NAME>-warnings.json`.
