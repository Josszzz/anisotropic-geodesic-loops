# geodesic-loops

Geodesic paths and closed geodesic loops on triangulated surface meshes, with
optional anisotropic (cardiac electrophysiology) metric.

Written in Rust, exposed as a Python extension via [PyO3](https://pyo3.rs) /
[maturin](https://www.maturin.rs).

## Algorithm

1. **Dijkstra** on the mesh graph initialises a shortest vertex-to-vertex path
   or non-contractible loop.
2. **Fan-unfolding shortening** (Polthier & Schmies 1998, Crane & Sharp
   flip-geodesics) iteratively straightens each interior vertex by unfolding the
   surrounding fan of triangles into the 2-D plane and replacing the vertex with
   the exact edge-crossing point(s) on the straight line between its neighbours.

Path points are represented as `(v0, v1, t)` triples:

| Type | v0 | v1 | t | 3-D position |
|------|----|----|---|--------------|
| Mesh vertex | `v` | `v` | `0.0` | `positions[v]` |
| Edge interior | `v0` | `v1` | `t ∈ (0,1)` | `positions[v0] + t·(positions[v1]−positions[v0])` |

## Install

```bash
cd geodesic_loops
pip install maturin
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop --release
```

Python ≥ 3.9, numpy required.

## Quick start

```python
import numpy as np
import geodesic_loops as gl

# Build a mesh
verts = np.array([...], dtype=np.float64)   # (N, 3)
faces = np.array([...], dtype=np.int64)     # (M, 3)
mesh  = gl.TriMesh(verts, faces)
# TriMesh(n_vertices=160, n_faces=280, n_edges=440)

# Point-to-point geodesic (Dijkstra + flip shortening)
path = gl.geodesic_path(mesh, src=0, dst=42)

polyline = path.polyline(mesh)   # (K, 3) float64 ndarray — 3-D coordinates
ep       = path.edge_params()    # list of (v0: int, v1: int, t: float)
length   = path.length(mesh)
stats    = path.stats(mesh)      # dict: length, n_vertices, is_closed, min/max angle

# Closed geodesic loop
loop = gl.geodesic_loop(mesh, seed=0)
loop.polyline(mesh)              # (K+1, 3) — first point repeated to close polygon
```

## API reference

### `TriMesh`

```python
mesh = gl.TriMesh(vertices, faces)   # vertices: (N,3) f64, faces: (M,3) i64

mesh.n_vertices   # int
mesh.n_faces      # int
mesh.n_edges      # int
mesh.edge_lengths()  # (n_edges,) f64 ndarray — Euclidean edge lengths
```

#### Metrics

```python
# Default: Euclidean
mesh.clear_metric()

# Isotropic conduction speed map
# Edge weight = |e| / harmonic_mean(c[v0], c[v1])
speeds = np.ones(mesh.n_vertices)
speeds[scar_verts] = 0.1           # 10× slower in scar region
mesh.set_isotropic_speeds(speeds)
# mesh.set_vertex_speeds(speeds)   # alias

# Anisotropic fiber-transverse metric (cardiac EP)
# g_v(d) = |d|²/c_t² + (d·f)²·(1/c_f² − 1/c_t²)
# Edge weight = sqrt( (g_v0(d) + g_v1(d)) / 2 )
fiber_dirs        = np.zeros((mesh.n_vertices, 3))
fiber_dirs[:, 2]  = 1.0                            # fibers along z-axis
speeds_fiber      = np.full(mesh.n_vertices, 1.0)  # fast along fiber
speeds_transverse = np.full(mesh.n_vertices, 0.3)  # slow transverse
mesh.set_fiber_metric(fiber_dirs, speeds_fiber, speeds_transverse)
```

### `GeodesicPath`

```python
path.n_vertices          # total number of path points (vertices + edge crossings)
path.is_closed           # bool

path.polyline(mesh)      # (K, 3) float64 ndarray; closed paths append first point
path.length(mesh)        # float — metric length
path.stats(mesh)         # dict with keys: length, n_vertices, is_closed,
                         #                 min_angle_rad, max_angle_rad

# Full edge-parameter representation (all points)
# Vertex(v)      → (v,  v,  0.0)
# Edge{v0,v1,t}  → (v0, v1, t)
path.edge_params()       # list of (v0: int, v1: int, t: float)

# Vertex-only indices (skips edge-interior points)
path.vertices()          # int64 ndarray
```

### Module functions

```python
# Dijkstra + flip shortening
gl.geodesic_path(mesh, src, dst,
                 max_iterations=10000, rel_tol=1e-8)   # → GeodesicPath | None

gl.geodesic_loop(mesh, seed,
                 max_iterations=10000, rel_tol=1e-8)   # → GeodesicPath | None

# Non-contractible loop in the homotopy class crossing a set of cut edges
gl.geodesic_loop_with_cut(mesh, seed, cut_edges,
                          max_iterations=10000, rel_tol=1e-8)  # → GeodesicPath | None

# Dijkstra only (no shortening)
gl.shortest_path(mesh, src, dst)        # → GeodesicPath | None
gl.shortest_loop(mesh, seed)            # → GeodesicPath | None

# Flip shortening on an existing path
gl.flip_shorten(mesh, path,
                max_iterations=10000, rel_tol=1e-8)    # → GeodesicPath
```

## Edge-crossing representation

On a coarse mesh the shortening algorithm replaces mesh-vertex waypoints with
exact edge-interior crossing points, producing a piecewise-linear path that is
straight within each triangle face.

```
Flat 6×6 grid  src=(0,0)  dst=(2,1)
Euclidean distance: 2.236068
Geodesic length:    2.236068   ← exact

edge_params:
  (  0,   0, 0.0000)  vertex          3D = [0.0,  0.0,  0.0]
  (  7,   6, 0.5000)  edge (7→6)      3D = [1.0,  0.5,  0.0]
  (  7,  12, 0.3333)  edge (7→12)     3D = [1.33, 0.67, 0.0]
  ( 13,  13, 0.0000)  vertex          3D = [2.0,  1.0,  0.0]
```

## Running the demo

```bash
python examples/demo.py
```
