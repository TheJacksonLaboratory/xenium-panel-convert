use std::str::FromStr;

use hdf5_metno::{
    Container, File, Group, H5Type,
    types::{FixedAscii, VarLenUnicode},
};
use ndarray::{Array1, ArrayView, Dimension};
use serde::Serialize;
use strum::VariantNames;

pub(super) fn read_container(file: &File, path: &str) -> Result<Container, ReadH5FieldError> {
    // If we can read the field as a group, great, if not, try reading as a dataset
    let container = match file.group(path) {
        Ok(g) => g.as_container().expect("a group should be a container"),
        Err(_) => file
            .dataset(path)
            .and_then(|ds| ds.as_container())
            .map_err(|err| {
                ReadH5FieldError::new_invalid_field(&err, file, path, FieldType::Container)
            })?,
    };

    Ok(container)
}

pub(super) fn read_attribute<T: H5Type>(
    container: &Container,
    path: &str,
) -> Result<T, ReadH5FieldError> {
    container
        .attr(path)
        .and_then(|a| a.read_scalar())
        .map_err(|err| {
            ReadH5FieldError::new_invalid_field(
                &err,
                &container.file().expect("file should be available"),
                path,
                FieldType::Attribute,
            )
        })
}

pub(super) fn read_dataset_raw<T: H5Type>(
    file: &File,
    path: &str,
) -> Result<Vec<T>, ReadH5FieldError> {
    file.dataset(path)
        .and_then(|ds| ds.read_raw())
        .map_err(|err| ReadH5FieldError::new_invalid_field(&err, file, path, FieldType::Dataset))
}

pub(super) fn read_1d_string_dataset(
    file: &File,
    path: &str,
) -> Result<Array1<VarLenUnicode>, ReadH5FieldError> {
    let encoding_type: VarLenUnicode =
        read_attribute(&read_container(file, path)?, "encoding-type")?;

    let encoding_type = StringEncodingType::from_str(&encoding_type).map_err(|_| {
        ReadH5FieldError::UnknownEncodingType {
            object_path: path.to_owned(),
            found: encoding_type.to_string(),
            expected: StringEncodingType::VARIANTS,
        }
    })?;

    match encoding_type {
        StringEncodingType::Categorical => read_categorical_array(file, path),
        StringEncodingType::StringArray => read_string_array(file, path),
        StringEncodingType::NullableStringArray => read_nullable_string_array(file, path),
    }
}

fn read_categorical_array(
    file: &File,
    path: &str,
) -> Result<Array1<VarLenUnicode>, ReadH5FieldError> {
    let mut null_indices = Vec::new();

    let codes = read_1d_dataset::<i32>(file, &format!("{path}/codes"))?;
    let categories = read_1d_dataset::<VarLenUnicode>(file, &format!("{path}/categories"))?;

    let array = codes
        .iter()
        .enumerate()
        .filter_map(|(i, code)| {
            if *code == -1 {
                null_indices.push(i);
                None
            } else {
                #[expect(clippy::cast_sign_loss)]
                Some(categories[*code as usize].clone())
            }
        })
        .collect();

    if null_indices.is_empty() {
        Ok(array)
    } else {
        Err(ReadH5FieldError::NullValues {
            indices: null_indices,
            object_path: path.to_owned(),
        })
    }
}

fn read_string_array(file: &File, path: &str) -> Result<Array1<VarLenUnicode>, ReadH5FieldError> {
    read_1d_dataset(file, path)
}

fn read_nullable_string_array(
    file: &File,
    path: &str,
) -> Result<Array1<VarLenUnicode>, ReadH5FieldError> {
    let is_null_array = read_1d_dataset::<bool>(file, &format!("{path}/mask"))?;
    let null_indices: Vec<_> = is_null_array
        .iter()
        .enumerate()
        .filter_map(|(i, is_null)| is_null.then_some(i))
        .collect();

    if null_indices.is_empty() {
        read_string_array(file, &format!("{path}/values"))
    } else {
        Err(ReadH5FieldError::NullValues {
            indices: null_indices,
            object_path: path.to_owned(),
        })
    }
}

fn read_1d_dataset<T: H5Type>(file: &File, path: &str) -> Result<Array1<T>, ReadH5FieldError> {
    file.dataset(path)
        .and_then(|ds| ds.read_1d())
        .map_err(|err| ReadH5FieldError::new_invalid_field(&err, file, path, FieldType::Dataset))
}

#[cfg(test)]
pub(crate) fn read_test_1d_dataset<T: H5Type>(
    file: &File,
    path: &str,
) -> Result<Array1<T>, ReadH5FieldError> {
    read_1d_dataset(file, path)
}

pub(super) fn to_ascii<const N: usize>(s: &VarLenUnicode) -> FixedAscii<N> {
    FixedAscii::from_ascii(&s).expect("all strings are ASCII in this context")
}

pub(super) fn create_h5_group(file: &File, path: &str) -> Result<Group, WriteH5ObjectError> {
    file.create_group(path).map_err(|e| WriteH5ObjectError {
        object_path: path.to_owned(),
        reason: e.to_string(),
    })
}

pub(super) fn write_dataset_to_h5_group<'d, A, T, D>(
    group: &Group,
    path: &str,
    data: A,
) -> Result<(), WriteH5ObjectError>
where
    A: Into<ArrayView<'d, T, D>>,
    T: H5Type,
    D: Dimension,
{
    group
        .new_dataset_builder()
        .with_data(data)
        .create(path)
        .map_err(|e| WriteH5ObjectError {
            object_path: path.to_owned(),
            reason: e.to_string(),
        })?;

    Ok(())
}

#[derive(Clone, Copy, Debug, strum::EnumString, strum::VariantNames)]
#[strum(serialize_all = "kebab-case")]
enum StringEncodingType {
    Categorical,
    StringArray,
    NullableStringArray,
}

#[derive(Debug, Clone, Serialize, thiserror::Error)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum ReadH5FieldError {
    #[error(
        "{object_path} could not be read as a {field_type} ({hdf5_error}) - ensure the correct \
         column name was provided (available objects in H5AD: {:?})",
        available_objects
    )]
    InvalidH5ObjectPath {
        hdf5_error: String,
        object_path: String,
        field_type: FieldType,
        available_objects: Vec<String>,
    },
    #[error(
        "null values were found at the given indices of {object_path} - ensure every element of \
         the array has a value"
    )]
    NullValues {
        indices: Vec<usize>,
        object_path: String,
    },
    #[error(
        "{object_path} has an unknown encoding type {found}, expected one of {expected:?} - \
         ensure the file was written by scanpy"
    )]
    UnknownEncodingType {
        object_path: String,
        found: String,
        expected: &'static [&'static str],
    },
}

impl ReadH5FieldError {
    pub(super) fn new_invalid_field(
        err: &hdf5_metno::Error,
        file: &File,
        object_path: &str,
        field_type: FieldType,
    ) -> Self {
        Self::InvalidH5ObjectPath {
            hdf5_error: err.to_string(),
            field_type,
            object_path: object_path.to_owned(),
            available_objects: recurse_through_group(file),
        }
    }
}

fn recurse_through_group(group: &Group) -> Vec<String> {
    group
        .iter_visit_default(
            Vec::with_capacity(32),
            |group: &Group, name, _, object_names| {
                if group.dataset(name).is_ok() {
                    object_names.push(format!("{}/{name}", group.name()));
                    return true;
                }

                if let Ok(nested_group) = group.group(name) {
                    let mut sub_object_names = recurse_through_group(&nested_group);
                    object_names.append(&mut sub_object_names);
                }

                true
            },
        )
        .expect("iterating over a file should not produce errors")
}

#[derive(Debug, Clone, Copy, Serialize, strum::Display)]
#[serde(rename_all = "snake_case")]
#[strum(serialize_all = "snake_case")]
pub enum FieldType {
    Attribute,
    Container,
    Dataset,
}

#[derive(Debug, Clone, Serialize)]
pub struct WriteH5ObjectError {
    pub object_path: String,
    pub reason: String,
}

#[cfg(test)]
mod tests {
    use hdf5_metno::{File, types::VarLenUnicode};
    use ndarray::Array1;

    use crate::reference_dataset::h5_util::{ReadH5FieldError, read_1d_string_dataset};

    fn read_obs_column(column: &str) -> Result<Array1<VarLenUnicode>, ReadH5FieldError> {
        let file = File::open("test-data/csr_adata.h5ad").unwrap();

        read_1d_string_dataset(&file, &format!("obs/{column}"))
    }

    #[test]
    fn unannotated_cells_are_rejected() {
        let ReadH5FieldError::NullValues {
            indices,
            object_path: _,
        } = read_obs_column("annotation_missing").unwrap_err()
        else {
            panic!("expecting null-values error");
        };

        assert_eq!(indices, [9]);
    }
}
