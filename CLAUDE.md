# CLAUDE.md — anisotropic-geodesic-loops

## Project overview

Compute **organizing cycles of cardiac reentry** as closed geodesics on a triangulated epicardial surface mesh.

The electrophysiological metric is `g_c = c⁻² g₀`, where `c(x)` is the local conduction velocity and `g₀` the Euclidean surface metric. The organizing cycle `γ` minimises `L_{g_c}(γ) = ∫_γ ds / c(x)`.

Implemented as a **Rust library** (`geodesic_loops/`) exposed to Python via [PyO3](https://pyo3.rs) / [maturin](https://www.maturin.rs).

---

## Build & install

```bash
# Pure Rust tests (no Python needed)
cargo test

# Python extension — dev/editable mode
cd geodesic_loops
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop --release

# Build wheel
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin build --features python
pip install target/wheels/geodesic_loops-*.whl --break-system-packages --force-reinstall
```

Requirements: Rust stable, Python ≥ 3.9, numpy, maturin.
`PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1` is required for Python ≥ 3.12.

---

## Codebase map

```
geodesic_loops/
  src/
    lib.rs              — module declarations, #[cfg(feature="python")] gating
    mesh.rs             — HalfedgeMesh, Metric enum (Euclidean/IsotropicSpeed/FiberTensor)
    geodesic.rs         — Dijkstra, shortest_loop_through, shortest_loop_crossing_cut
    flip_geodesics.rs   — PathPoint, GeodesicPath, fan-unfolding shortening (shorten_path)
    python_bindings.rs  — PyO3 bindings: TriMesh, GeodesicPath, module functions
    tests.rs            — 25 unit tests
  Cargo.toml            — pyo3/numpy optional under [features] python
  pyproject.toml        — maturin config
examples/
  demo.py               — end-to-end Python demo
paper_draft.md          — theory: cardiac reentry, eikonal equation, Riemannian metric
implementation.md       — optimisation formulation (FR): variational energy J(u,γ)
```

---

## Algorithms

| Step | Location | Description |
|------|----------|-------------|
| Dijkstra | `geodesic.rs` | Shortest path on mesh graph; edge weights from active metric |
| Fan unfolding | `flip_geodesics.rs` | Polthier & Schmies 1998 / Crane & Sharp: replaces Vertex waypoints with exact edge-crossing points |
| Non-contractible loops | `geodesic.rs` | Doubled-graph (vertex × parity) via `shortest_loop_crossing_cut` |

---

## Key design decisions & gotchas

- **Python feature flag**: `[features] python = ["dep:pyo3", ...]` — pure Rust tests never link libpython.
- **PathPoint** = `Vertex(usize)` | `Edge { v0, v1, t: f64 }`. Closed loops do **not** repeat the first vertex; the closing segment is implicit.
- **`metric_edge_len_between`** returns `f64::INFINITY` for non-adjacent vertices — never fall back to Euclidean distance.
- **Fan unfolding**: the geodesic crosses **radial edges** (B → fan vertex at origin), NOT the outer fan edges (V_k → V_{k+1}). Multiple crossings per fan are possible.
- **`flat_grid` test mesh**: alternating diagonals (even i+j → a-d diagonal; odd → b-c) for near-isotropic behaviour.
- **`shorten_path` loop**: when `straighten_at_vertex` returns N replacements, splice them in-place via `pts.splice(cur_i..=cur_i, replacements)` and do not advance `i`.

---

## Python API

```python
import geodesic_loops as gl

# Mesh
mesh = gl.TriMesh(verts, faces)               # verts: (N,3) f64, faces: (M,3) i64
mesh.set_isotropic_speeds(speeds)             # (N,) f64 — travel-time metric
mesh.set_fiber_metric(dirs, c_f, c_t)         # (N,3), (N,), (N,) — anisotropic EP metric
mesh.clear_metric()                           # reset to Euclidean

# Paths
path = gl.geodesic_path(mesh, src, dst)       # Dijkstra + flip shortening → GeodesicPath | None
loop = gl.geodesic_loop(mesh, seed)           # closed geodesic loop → GeodesicPath | None
loop = gl.geodesic_loop_with_cut(mesh, seed, cut_edges)  # non-contractible loop

# GeodesicPath methods
path.polyline(mesh)    # (K, 3) ndarray — 3-D coords of all path points
path.edge_params()     # list of (v0, v1, t): Vertex(v) → (v,v,0.0), Edge → (v0,v1,t∈(0,1))
path.vertices()        # int64 ndarray — Vertex-only indices
path.length(mesh)      # float — metric length
path.stats(mesh)       # dict: length, n_vertices, is_closed, min_angle_rad, max_angle_rad
path.n_vertices        # total point count (vertices + edge crossings)
path.is_closed         # bool
```

---

## Tests

```bash
cargo test              # 25 tests, all must pass
python examples/demo.py # end-to-end smoke test
```

Key regression: `flip_tests::shortening_produces_edge_crossings` — path `[0,7,13]` on 6×6 flat grid must produce ≥1 `Edge` crossing and converge to length ≈ √5 ≈ 2.2361.

---

## Scientific context (broader project)

Joint reconstruction from observed activation maps `τ_obs`:

- **conduction field** `c(x)` (parameterised as `u = log c`)
- **organizing cycle** `γ` (closed geodesic in homotopy class α)

Energy: `J(u,γ) = J_data + λ J_reg + μ J_geo + ν J_phys`, minimised by alternating optimisation. See `implementation.md` (methods, FR) and `paper_draft.md` (theory).
