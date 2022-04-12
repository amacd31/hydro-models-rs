use numpy::{IntoPyArray, PyArray1, PyReadonlyArrayDyn};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::models::gr4j::{GR4JModel, GR4JParams};

#[pymodule]
pub fn hydromodels(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<GR4JModel>()?;
    m.add_class::<GR4JParams>()?;
    Ok(())
}

#[pymethods]
impl GR4JParams {
    #[new]
    fn py_new(x1: f64, x2: f64, x3: f64, x4: f64) -> PyResult<Self> {
        match GR4JParams::new(x1, x2, x3, x4) {
            Ok(params) => Ok(params),
            Err(e) => Err(PyValueError::new_err(e)),
        }
    }
}

#[pymethods]
impl GR4JModel {
    #[new]
    fn new(params: GR4JParams) -> Self {
        GR4JModel::create(params)
    }

    fn run_model<'py>(
        &mut self,
        py: Python<'py>,
        precip: PyReadonlyArrayDyn<f64>,
        pet: PyReadonlyArrayDyn<f64>,
    ) -> &'py PyArray1<f64> {
        let precip = precip.as_slice().unwrap();
        let pet = pet.as_slice().unwrap();

        self.run(&precip, &pet).into_pyarray(py)
    }
}
