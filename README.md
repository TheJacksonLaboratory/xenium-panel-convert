# xenium-panel-convert

A command-line utility for converting files to the formats accepted by the 10x Genomics Xenium Panel Designer.

This command-line tool validates and converts two types of files to formats the [10x Genomics Xenium Panel Designer](https://www.10xgenomics.com/support/software/xenium-panel-designer/latest) accepts. These file types are:

- CSV-formatted gene-lists
- [scanpy](https://github.com/scverse/scanpy)-generated anndata files (H5AD)

Run `xp-convert --help` for more information.

## Installation

Install the latest version from the [releases page](https://github.com/TheJacksonLaboratory/xenium-panel-convert/releases). The installer script will also install a sibling binary called `xenium-panel-convert-update` - run that to upgrade to the latest version.
