use std::{fs, io};

use anyhow::{Context, ensure};
use camino::Utf8Path;
use serde::Serialize;

pub(crate) fn write_json_to_file(data: &impl Serialize, path: &Utf8Path) -> anyhow::Result<()> {
    ensure!(!path.exists(), "cannot overwrite file at {path}");

    let error_message = || format!("failed to write error report to {path}");
    let file = fs::File::create(path).with_context(error_message)?;

    serde_json::to_writer(io::BufWriter::new(file), data).with_context(error_message)?;

    Ok(())
}

pub(crate) fn write_csv_with_header_to_file(
    header: &[&str],
    rows: &[impl Serialize],
    path: &Utf8Path,
) -> anyhow::Result<()> {
    ensure!(!path.exists(), "cannot overwrite file at {path}");

    let mut writer = csv::WriterBuilder::new().has_headers(false).from_path(path)?;

    writer.write_record(header)?;

    for row in rows {
        writer.serialize(row)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use camino::Utf8Path;
    use tempfile::TempDir;

    use crate::write::write_json_to_file;

    #[test]
    fn does_not_overwrite() {
        let dir = TempDir::new().unwrap();
        let path = Utf8Path::from_path(dir.path()).unwrap().join("output.json");

        write_json_to_file(&"first", &path).unwrap();

        assert!(
            write_json_to_file(&"second", &path).is_err(),
            "an existing file should not be overwritten"
        );
    }
}
