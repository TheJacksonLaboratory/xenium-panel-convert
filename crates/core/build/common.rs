use std::{
    fs,
    io::{self, Write},
    path::Path,
};

use anyhow::Context;
use serde::de::DeserializeOwned;

pub(crate) fn parse_gene_list_from_csv<T: DeserializeOwned>(raw_gene_list: &[u8]) -> anyhow::Result<Vec<T>> {
    let mut reader = csv::ReaderBuilder::new().comment(Some(b'#')).from_reader(raw_gene_list);

    Ok(reader.deserialize().collect::<Result<_, _>>()?)
}

pub(crate) fn write_map_to_file(path: &Path, map_name: &str, map: &phf_codegen::Map<'_, &str>) -> anyhow::Result<()> {
    let file = fs::File::create(path).with_context(|| format!("failed to write file {}", path.to_str().unwrap()))?;
    let mut file_writer = io::BufWriter::new(file);

    writeln!(
        file_writer,
        "pub(super) static {map_name}: phf::Map<&'static str, &'static str> = {};",
        map.build()
    )?;

    Ok(())
}
