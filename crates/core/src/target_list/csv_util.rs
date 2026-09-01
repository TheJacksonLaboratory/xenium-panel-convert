use std::collections::HashMap;

use csv::StringRecord;

pub(super) fn read_csv_trimmed(target_list: &str) -> csv::Reader<&[u8]> {
    let target_list = target_list.trim();

    csv::ReaderBuilder::new()
        .trim(csv::Trim::All)
        .from_reader(target_list.as_bytes())
}

pub(super) fn rename_fields(
    original_fieldnames: &StringRecord,
    field_aliases: &HashMap<&str, &str>,
) -> StringRecord {
    let mut renamed_fields = StringRecord::new();

    for original in original_fieldnames {
        let renamed = field_aliases.get(original).unwrap_or(&original);

        renamed_fields.push_field(renamed);
    }

    renamed_fields
}

#[cfg(test)]
mod tests {
    use crate::target_list::csv_util::rename_fields;

    #[test]
    fn renaming_fields() {
        let original_fieldnames = ["field1", "field2"].iter().collect();
        let field_aliases = [("field1", "field_1")].into_iter().collect();

        let renamed_fields = rename_fields(&original_fieldnames, &field_aliases);

        assert_eq!(
            renamed_fields,
            ["field_1", "field2"][..],
            "failed to rename fields"
        );
    }
}
