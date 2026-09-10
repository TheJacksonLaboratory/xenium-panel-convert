use camino::Utf8Path;
use xenium_panel_convert_core::all::validate_target_list_and_reference_dataset_compatibility;

use crate::{
    TargetListCliOptions, convert_target_list, dataset_name,
    reference_datasets::{ReferenceDatasetCliOptions, convert_reference_dataset_and_write},
    write::write_json_to_file,
};

pub(super) fn convert_target_list_and_reference_datasets(
    targets_options: &TargetListCliOptions,
    references_options: &ReferenceDatasetCliOptions,
    output_dir: &Utf8Path,
) -> anyhow::Result<()> {
    let target_list = convert_target_list(targets_options, output_dir)?;

    let Some(target_list) = target_list else {
        return Ok(());
    };

    for spec in &references_options.reference_datasets {
        let Ok(ds) = convert_reference_dataset_and_write(spec, output_dir) else {
            continue;
        };

        let warnings = validate_target_list_and_reference_dataset_compatibility(
            target_list.targets(),
            targets_options.species,
            &ds,
            spec.transcriptome_name,
            spec.flex,
        );

        if warnings.is_empty() {
            continue;
        }

        let ds_name = dataset_name(&spec.path, spec.rename.as_deref())?;
        let warnings_path = format!("{ds_name}-warnings.json");
        write_json_to_file(&warnings, &output_dir.join(warnings_path))?;
    }

    Ok(())
}
