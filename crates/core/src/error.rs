use std::fmt::Display;

use serde::Serialize;

/// An error paired with the message written to the JSON error reports.
///
/// `hint` is the inner error's [`Display`] output. Every `#[error(...)]`
/// message in this crate is user-facing, and must describe what went wrong
/// before giving the remedy after a dash - for example, "the field 'group' is
/// missing - add it to the CSV".
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct Hinted<E> {
    #[serde(flatten)]
    pub error: E,
    pub hint: String,
}

impl<E: Display> Hinted<E> {
    pub(crate) fn new(error: E) -> Self {
        Self {
            hint: error.to_string(),
            error,
        }
    }
}

pub(crate) fn collect_error<T, E1, E2>(result: Result<T, E1>, errors: &mut Vec<E2>) -> Option<T>
where
    E1: Into<E2>,
{
    match result {
        Ok(value) => Some(value),
        Err(error) => {
            errors.push(error.into());

            None
        }
    }
}
