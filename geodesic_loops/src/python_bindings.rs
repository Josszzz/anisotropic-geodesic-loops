use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};

use crate::mesh::HalfedgeMesh;
use crate::flip_geodesics::{FlipConfig, GeodesicPath};

// ---- Python wrapper for HalfedgeMesh ----

/// A triangulated surface mesh.
#[pyclass(name = "TriMesh")]
pub struct PyTriMesh {
    pub(crate) mesh: HalfedgeMesh,
}

#[pymethods]
impl PyTriMesh {
    #[new]
    fn new(
        _py: Python<'_>,
        vertices: PyReadonlyArray2<f64>,
        faces: PyReadonlyArray2<i64>,
    ) -> PyResult<Self> {
        let verts_arr = vertices.as_array();
        let faces_arr = faces.as_array();

        if verts_arr.ncols() != 3 {
            return Err(PyValueError::new_err("vertices must have shape (N, 3)"));
        }
        if faces_arr.ncols() != 3 {
            return Err(PyValueError::new_err("faces must have shape (M, 3)"));
        }

        let n_v = verts_arr.nrows();
        let n_f = faces_arr.nrows();

        let verts: Vec<[f64; 3]> = (0..n_v)
            .map(|i| [verts_arr[[i, 0]], verts_arr[[i, 1]], verts_arr[[i, 2]]])
            .collect();
        let tris: Vec<[usize; 3]> = (0..n_f)
            .map(|i| {
                [faces_arr[[i, 0]] as usize,
                 faces_arr[[i, 1]] as usize,
                 faces_arr[[i, 2]] as usize]
            })
            .collect();

        Ok(Self { mesh: HalfedgeMesh::from_vertices_faces(&verts, &tris) })
    }

    /// Reset to Euclidean metric (default).
    fn clear_metric(&mut self) { self.mesh.clear_metric(); }

    /// Set isotropic per-vertex conduction speed map.
    ///
    /// Edge weight: ``|e| / harmonic_mean(c[v0], c[v1])`` (travel time).
    fn set_isotropic_speeds(&mut self, speeds: PyReadonlyArray1<f64>) -> PyResult<()> {
        let s = speeds.as_array();
        if s.len() != self.mesh.n_verts() {
            return Err(PyValueError::new_err(format!(
                "speeds length {} != n_verts {}", s.len(), self.mesh.n_verts()
            )));
        }
        self.mesh.set_isotropic_speeds(s.to_vec());
        Ok(())
    }

    /// Alias for ``set_isotropic_speeds`` (backward compatibility).
    fn set_vertex_speeds(&mut self, speeds: PyReadonlyArray1<f64>) -> PyResult<()> {
        self.set_isotropic_speeds(speeds)
    }

    /// Alias for ``clear_metric`` (backward compatibility).
    fn clear_vertex_speeds(&mut self) { self.mesh.clear_metric(); }

    /// Set anisotropic fiber-transverse metric for cardiac electrophysiology.
    ///
    /// Parameters
    /// ----------
    /// fiber_dirs : ndarray, shape (N, 3)
    ///     Per-vertex fiber direction unit vectors (tangent to the surface).
    /// speeds_fiber : ndarray, shape (N,)
    ///     Conduction speed along the fiber direction at each vertex.
    /// speeds_transverse : ndarray, shape (N,)
    ///     Conduction speed transverse to the fiber direction at each vertex.
    ///
    /// The metric quadratic form at vertex v is:
    ///   g_v(d) = |d|²/c_t(v)² + (d·f(v))² · (1/c_f(v)² − 1/c_t(v)²)
    ///
    /// Edge weight: sqrt( (g_v0(d) + g_v1(d)) / 2 )
    fn set_fiber_metric(
        &mut self,
        fiber_dirs: PyReadonlyArray2<f64>,
        speeds_fiber: PyReadonlyArray1<f64>,
        speeds_transverse: PyReadonlyArray1<f64>,
    ) -> PyResult<()> {
        let n = self.mesh.n_verts();
        let fd = fiber_dirs.as_array();
        let sf = speeds_fiber.as_array();
        let st = speeds_transverse.as_array();

        if fd.nrows() != n || fd.ncols() != 3 {
            return Err(PyValueError::new_err(format!(
                "fiber_dirs must have shape ({}, 3), got ({}, {})", n, fd.nrows(), fd.ncols()
            )));
        }
        if sf.len() != n {
            return Err(PyValueError::new_err(format!(
                "speeds_fiber length {} != n_verts {}", sf.len(), n
            )));
        }
        if st.len() != n {
            return Err(PyValueError::new_err(format!(
                "speeds_transverse length {} != n_verts {}", st.len(), n
            )));
        }

        let fiber_dirs_vec: Vec<[f64; 3]> = (0..n)
            .map(|i| [fd[[i, 0]], fd[[i, 1]], fd[[i, 2]]])
            .collect();
        self.mesh.set_fiber_metric(fiber_dirs_vec, sf.to_vec(), st.to_vec());
        Ok(())
    }

    #[getter] fn n_vertices(&self) -> usize { self.mesh.n_verts() }
    #[getter] fn n_faces(&self) -> usize { self.mesh.n_faces() }
    #[getter] fn n_edges(&self) -> usize { self.mesh.n_edges() }

    fn edge_lengths<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec(py, self.mesh.edge_lengths.clone())
    }

    fn __repr__(&self) -> String {
        format!("TriMesh(n_vertices={}, n_faces={}, n_edges={})",
            self.mesh.n_verts(), self.mesh.n_faces(), self.mesh.n_edges())
    }
}

// ---- Python wrapper for GeodesicPath ----

#[pyclass(name = "GeodesicPath")]
pub struct PyGeodesicPath {
    pub(crate) path: GeodesicPath,
}

#[pymethods]
impl PyGeodesicPath {
    /// Vertex indices for pure-Vertex path points as int64 array.
    /// Edge-crossing points are skipped; use `edge_params()` for the full path.
    fn vertices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i64>> {
        let v: Vec<i64> = self.path.vertex_indices().iter().map(|&x| x as i64).collect();
        PyArray1::from_vec(py, v)
    }

    /// Edge-parameter representation of all path points as an (N, 3) int/float array.
    ///
    /// Each row is ``[v0, v1, t]`` where the point is at
    /// ``positions[v0] + t * (positions[v1] - positions[v0])``.
    /// For a Vertex point, ``v0 == v1`` and ``t == 0``.
    ///
    /// Returns a Python list of (v0: int, v1: int, t: float) tuples.
    fn edge_params<'py>(&self, py: Python<'py>) -> PyObject {
        let params = self.path.edge_params();
        let list = pyo3::types::PyList::empty(py);
        for (v0, v1, t) in params {
            let tup = pyo3::types::PyTuple::new(py, &[
                v0.into_pyobject(py).unwrap().into_any().unbind(),
                v1.into_pyobject(py).unwrap().into_any().unbind(),
                t.into_pyobject(py).unwrap().into_any().unbind(),
            ]).unwrap();
            list.append(tup).unwrap();
        }
        list.into()
    }

    #[getter] fn is_closed(&self) -> bool { self.path.is_closed }
    #[getter] fn n_vertices(&self) -> usize { self.path.len() }

    /// 3-D polyline as ndarray of shape (K, 3).
    fn polyline<'py>(&self, py: Python<'py>, mesh: &PyTriMesh) -> Bound<'py, PyArray2<f64>> {
        let pts = self.path.to_polyline(&mesh.mesh);
        let data: Vec<Vec<f64>> = pts.iter().map(|p| p.to_vec()).collect();
        PyArray2::from_vec2(py, &data).unwrap()
    }

    /// Metric length of the path.
    fn length(&self, mesh: &PyTriMesh) -> f64 { self.path.metric_length(&mesh.mesh) }

    /// Statistics dict: length, n_vertices, is_closed, min_angle, max_angle (radians).
    fn stats<'py>(&self, py: Python<'py>, mesh: &PyTriMesh) -> PyObject {
        let s = crate::flip_geodesics::compute_stats(&mesh.mesh, &self.path);
        let d = pyo3::types::PyDict::new(py);
        d.set_item("length", s.length).unwrap();
        d.set_item("n_vertices", s.n_points).unwrap();
        d.set_item("is_closed", s.is_closed).unwrap();
        d.set_item("min_angle_rad", s.min_angle).unwrap();
        d.set_item("max_angle_rad", s.max_angle).unwrap();
        d.into()
    }

    fn __repr__(&self) -> String {
        format!("GeodesicPath(n_vertices={}, is_closed={})", self.path.len(), self.path.is_closed)
    }
}

// ---- module-level functions ----

/// Shortest point-to-point path (Dijkstra, no shortening).
#[pyfunction]
fn shortest_path(mesh: &PyTriMesh, src: usize, dst: usize) -> Option<PyGeodesicPath> {
    crate::geodesic::shortest_path(&mesh.mesh, src, dst)
        .map(|v| PyGeodesicPath { path: GeodesicPath::open(v) })
}

/// Shortest non-contractible loop through a seed vertex (Dijkstra only).
#[pyfunction]
fn shortest_loop(mesh: &PyTriMesh, seed: usize) -> Option<PyGeodesicPath> {
    crate::geodesic::shortest_loop_through(&mesh.mesh, seed)
        .map(|v| PyGeodesicPath { path: GeodesicPath::closed(v) })
}

/// Shorten an existing path/loop with the flip algorithm.
#[pyfunction]
#[pyo3(signature = (mesh, path, max_iterations=10000, rel_tol=1e-8))]
fn flip_shorten(
    mesh: &PyTriMesh, path: &PyGeodesicPath,
    max_iterations: usize, rel_tol: f64,
) -> PyGeodesicPath {
    let config = FlipConfig { max_iterations, rel_tol, angle_eps: 1e-6 };
    PyGeodesicPath { path: crate::flip_geodesics::shorten_path(&mesh.mesh, &path.path, &config) }
}

/// Dijkstra + flip shortening: point-to-point geodesic.
#[pyfunction]
#[pyo3(signature = (mesh, src, dst, max_iterations=10000, rel_tol=1e-8))]
fn geodesic_path(
    mesh: &PyTriMesh, src: usize, dst: usize,
    max_iterations: usize, rel_tol: f64,
) -> Option<PyGeodesicPath> {
    let config = FlipConfig { max_iterations, rel_tol, angle_eps: 1e-6 };
    crate::flip_geodesics::geodesic_path(&mesh.mesh, src, dst, &config)
        .map(|p| PyGeodesicPath { path: p })
}

/// Dijkstra + flip shortening: closed geodesic loop.
#[pyfunction]
#[pyo3(signature = (mesh, seed, max_iterations=10000, rel_tol=1e-8))]
fn geodesic_loop(
    mesh: &PyTriMesh, seed: usize,
    max_iterations: usize, rel_tol: f64,
) -> Option<PyGeodesicPath> {
    let config = FlipConfig { max_iterations, rel_tol, angle_eps: 1e-6 };
    crate::flip_geodesics::geodesic_loop(&mesh.mesh, seed, &config)
        .map(|p| PyGeodesicPath { path: p })
}

/// Geodesic loop in a homotopy class defined by a list of cut edge indices.
#[pyfunction]
#[pyo3(signature = (mesh, seed, cut_edges, max_iterations=10000, rel_tol=1e-8))]
fn geodesic_loop_with_cut(
    mesh: &PyTriMesh, seed: usize, cut_edges: Vec<usize>,
    max_iterations: usize, rel_tol: f64,
) -> Option<PyGeodesicPath> {
    let cut_set: std::collections::HashSet<usize> = cut_edges.into_iter().collect();
    let config = FlipConfig { max_iterations, rel_tol, angle_eps: 1e-6 };
    crate::flip_geodesics::geodesic_loop_with_cut(&mesh.mesh, seed, &cut_set, &config)
        .map(|p| PyGeodesicPath { path: p })
}

// ---- module definition ----

/// Geodesic loops on triangulated meshes.
///
/// Computes point-to-point geodesics and closed geodesic loops on triangle
/// meshes, with optional anisotropic (electrophysiological) metric. Paths are
/// initialized by Dijkstra and shortened by the Crane & Sharp flip algorithm.
///
/// Quick start::
///
///     import numpy as np
///     import geodesic_loops as gl
///
///     # Build a mesh
///     verts = np.array([...], dtype=np.float64)  # (N, 3)
///     faces = np.array([...], dtype=np.int64)     # (M, 3)
///     mesh = gl.TriMesh(verts, faces)
///
///     # Point-to-point geodesic
///     path = gl.geodesic_path(mesh, src=0, dst=42)
///     pts = path.polyline(mesh)  # (K, 3) ndarray
///
///     # Closed geodesic loop
///     loop = gl.geodesic_loop(mesh, seed=0)
///
///     # With anisotropic metric
///     speeds = np.ones(mesh.n_vertices)
///     speeds[scar_region] = 0.1  # slow zone
///     mesh.set_vertex_speeds(speeds)
///     loop_aniso = gl.geodesic_loop(mesh, seed=0)
#[pymodule]
pub fn geodesic_loops(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTriMesh>()?;
    m.add_class::<PyGeodesicPath>()?;
    m.add_function(wrap_pyfunction!(shortest_path, m)?)?;
    m.add_function(wrap_pyfunction!(shortest_loop, m)?)?;
    m.add_function(wrap_pyfunction!(flip_shorten, m)?)?;
    m.add_function(wrap_pyfunction!(geodesic_path, m)?)?;
    m.add_function(wrap_pyfunction!(geodesic_loop, m)?)?;
    m.add_function(wrap_pyfunction!(geodesic_loop_with_cut, m)?)?;
    Ok(())
}
