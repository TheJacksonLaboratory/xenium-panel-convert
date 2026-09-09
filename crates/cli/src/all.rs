use camino::Utf8Path;
use xenium_panel_convert_core::all::validate_target_list_and_reference_dataset_compatibility;

use crate::{
    TargetListCliOptions, convert_reference_datasets, convert_target_list, dataset_name,
    reference_datasets::ReferenceDatasetCliOptions, write::write_json_to_file,
};

pub(super) fn convert_target_list_and_reference_datasets(
    targets_options: &TargetListCliOptions,
    references_options: &ReferenceDatasetCliOptions,
    output_dir: &Utf8Path,
) -> anyhow::Result<()> {
    let target_list = convert_target_list(targets_options, output_dir)?;
    let reference_datasets = convert_reference_datasets(references_options, output_dir)?;

    let Some(target_list) = target_list else {
        return Ok(());
    };

    for (ds, options) in &reference_datasets {
        let warnings = validate_target_list_and_reference_dataset_compatibility(
            target_list.targets(),
            targets_options.species,
            ds,
            options.transcriptome_name,
            options.flex,
        );

        if warnings.is_empty() {
            continue;
        }

        let ds_name = dataset_name(&options.path, options.rename.as_deref())?;
        let warnings_path = format!("{ds_name}-warnings.json");
        write_json_to_file(&warnings, &output_dir.join(warnings_path))?;
    }

    Ok(())
}
