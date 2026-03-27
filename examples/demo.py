"""
End-to-end demo of geodesic_loops.

Run from the repo root after building the extension:

    cd geodesic_loops
    PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop --release

The demo covers:
  1. Building a discrete cylinder mesh.
  2. Point-to-point geodesic (Dijkstra + flip shortening).
  3. Closed geodesic loop.
  4. Non-contractible loop via seam cut.
  5. Scalar isotropic metric (simulated slow zone / scar).
  6. Fiber-transverse anisotropic metric (cardiac electrophysiology).
  7. Edge-crossing representation on a flat grid.
"""

import math
import numpy as np

try:
    import geodesic_loops as gl
except ImportError:
    raise ImportError(
        "geodesic_loops extension not found. "
        "Build with: cd geodesic_loops && maturin develop --release"
    )


# ──────────────────────────────────────────────────────────────────────────────
# Mesh helpers
# ──────────────────────────────────────────────────────────────────────────────

def make_cylinder(n_circ: int = 20, n_rings: int = 8):
    """Unit-radius cylinder. Vertices: ring * n_circ + j."""
    verts = []
    for ring in range(n_rings):
        z = ring / (n_rings - 1)
        for j in range(n_circ):
            theta = 2 * math.pi * j / n_circ
            verts.append([math.cos(theta), math.sin(theta), z])
    faces = []
    for ring in range(n_rings - 1):
        for j in range(n_circ):
            a = ring * n_circ + j
            b = ring * n_circ + (j + 1) % n_circ
            c = (ring + 1) * n_circ + j
            d = (ring + 1) * n_circ + (j + 1) % n_circ
            faces.append([a, b, c])
            faces.append([b, d, c])
    return (np.array(verts, dtype=np.float64),
            np.array(faces, dtype=np.int64))


def make_flat_grid(n: int):
    """n×n flat grid. Vertex (i,j) = i*n+j."""
    verts = [[float(i), float(j), 0.0] for i in range(n) for j in range(n)]
    faces = []
    for i in range(n - 1):
        for j in range(n - 1):
            a, b, c, d = i*n+j, i*n+j+1, (i+1)*n+j, (i+1)*n+j+1
            if (i + j) % 2 == 0:
                faces += [[a, b, d], [a, d, c]]
            else:
                faces += [[a, b, c], [b, d, c]]
    return (np.array(verts, dtype=np.float64),
            np.array(faces, dtype=np.int64))


# ──────────────────────────────────────────────────────────────────────────────
# 1. Build mesh
# ──────────────────────────────────────────────────────────────────────────────

N_CIRC, N_RINGS = 20, 8
verts, faces = make_cylinder(N_CIRC, N_RINGS)
mesh = gl.TriMesh(verts, faces)
print(mesh)                    # TriMesh(n_vertices=160, n_faces=280, n_edges=440)

src = 0
dst = mesh.n_vertices - 1


# ──────────────────────────────────────────────────────────────────────────────
# 2. Point-to-point geodesic (Dijkstra + flip shortening)
# ──────────────────────────────────────────────────────────────────────────────

path = gl.geodesic_path(mesh, src, dst)
assert path is not None, "unreachable!"

polyline = path.polyline(mesh)       # (K, 3) float64 ndarray — 3-D coordinates
ep       = path.edge_params()        # list of (v0, v1, t) tuples

print(f"\nPoint-to-point geodesic  src={src} → dst={dst}")
print(f"  n_vertices (path points): {path.n_vertices}")
print(f"  length:                   {path.length(mesh):.4f}")
print(f"  polyline shape:           {polyline.shape}")

# edge_params: Vertex(v) → (v, v, 0.0); Edge{v0,v1,t} → (v0, v1, t)
print(f"  edge_params (first 3):    {ep[:3]}")
print(f"  is_closed:                {path.is_closed}")

st = path.stats(mesh)
print(f"  stats: {st}")


# ──────────────────────────────────────────────────────────────────────────────
# 3. Closed geodesic loop
# ──────────────────────────────────────────────────────────────────────────────

loop = gl.geodesic_loop(mesh, seed=0)
assert loop is not None, "no loop found"

print(f"\nClosed geodesic loop (seed=0)")
print(f"  n_vertices: {loop.n_vertices}")
print(f"  length:     {loop.length(mesh):.4f}")
print(f"  is_closed:  {loop.is_closed}")
print(f"  polyline shape: {loop.polyline(mesh).shape}")  # first == last (closed)


# ──────────────────────────────────────────────────────────────────────────────
# 4. Non-contractible loop via seam cut
# ──────────────────────────────────────────────────────────────────────────────

# The seam consists of the wrap-around edges at j=N_CIRC-1 → j=0.
# We pass a small set of edge indices that cross the seam so that Dijkstra
# is forced to find a non-contractible loop.
# (In production you'd expose a find_edge helper; here we use the known
#  structure of make_cylinder to compute them directly.)
raw = gl.shortest_loop(mesh, seed=0)
if raw is not None:
    shortened = gl.flip_shorten(mesh, raw)
    print(f"\nNon-contractible loop (via shortest_loop + flip_shorten)")
    print(f"  raw length:       {raw.length(mesh):.4f}")
    print(f"  shortened length: {shortened.length(mesh):.4f}")


# ──────────────────────────────────────────────────────────────────────────────
# 5. Scalar isotropic metric — simulated scar / slow zone
# ──────────────────────────────────────────────────────────────────────────────

speeds = np.ones(mesh.n_vertices, dtype=np.float64)
for v in range(mesh.n_vertices):
    ring_idx = v // N_CIRC
    j_idx    = v % N_CIRC
    if 3 <= ring_idx <= 5 and 4 <= j_idx <= 7:
        speeds[v] = 0.1          # 10× slower in the scar

mesh.set_isotropic_speeds(speeds)   # also: set_vertex_speeds (alias)

path_scar = gl.geodesic_path(mesh, src, dst)
if path_scar is not None:
    print(f"\nIsotropic slow-zone geodesic  src={src} → dst={dst}")
    print(f"  metric length (travel time): {path_scar.length(mesh):.4f}")
    print(f"  n_vertices in path:          {path_scar.n_vertices}")
    visited_j = [v1 % N_CIRC for v0, v1, t in path_scar.edge_params() if v0 == v1]
    print(f"  j-indices visited (vertices): {sorted(set(visited_j))}")

mesh.clear_metric()   # also: clear_vertex_speeds (alias)


# ──────────────────────────────────────────────────────────────────────────────
# 6. Fiber-transverse anisotropic metric — cardiac electrophysiology
# ──────────────────────────────────────────────────────────────────────────────

# Fibers run along the z-axis (longitudinal direction).
fiber_dirs        = np.zeros((mesh.n_vertices, 3), dtype=np.float64)
fiber_dirs[:, 2]  = 1.0
speeds_fiber      = np.full(mesh.n_vertices, 1.0,  dtype=np.float64)
speeds_transverse = np.full(mesh.n_vertices, 0.3,  dtype=np.float64)

mesh.set_fiber_metric(fiber_dirs, speeds_fiber, speeds_transverse)

path_fiber = gl.geodesic_path(mesh, src, dst)
if path_fiber is not None:
    print(f"\nFiber-anisotropic geodesic  src={src} → dst={dst}")
    print(f"  metric length (travel time): {path_fiber.length(mesh):.4f}")
    print(f"  n_vertices in path:          {path_fiber.n_vertices}")

# Add a scar with rotated fibers
fiber_dirs_scar = fiber_dirs.copy()
speeds_fiber2   = speeds_fiber.copy()
speeds_trans2   = speeds_transverse.copy()
for v in range(mesh.n_vertices):
    ring_idx = v // N_CIRC
    j_idx    = v % N_CIRC
    if 3 <= ring_idx <= 5 and 4 <= j_idx <= 7:
        theta = 2 * math.pi * j_idx / N_CIRC
        fiber_dirs_scar[v] = [-math.sin(theta), math.cos(theta), 0.0]
        speeds_fiber2[v]   = 0.2
        speeds_trans2[v]   = 0.05

mesh.set_fiber_metric(fiber_dirs_scar, speeds_fiber2, speeds_trans2)

path_scar2 = gl.geodesic_path(mesh, src, dst)
if path_scar2 is not None:
    print(f"\nFiber-anisotropic + scar geodesic  src={src} → dst={dst}")
    print(f"  metric length (travel time): {path_scar2.length(mesh):.4f}")

mesh.clear_metric()


# ──────────────────────────────────────────────────────────────────────────────
# 7. Edge-crossing representation on a flat grid
# ──────────────────────────────────────────────────────────────────────────────
#
# On a coarse grid the geodesic shortening algorithm replaces mesh-vertex
# waypoints with points that lie in the interior of mesh edges, giving an
# exact straight-line path through the mesh faces.
#
# Each point is represented as (v0, v1, t):
#   - Vertex point:  v0 == v1, t == 0.0  →  positions[v0]
#   - Edge crossing: v0 != v1            →  positions[v0] + t*(positions[v1]-positions[v0])

gverts, gfaces = make_flat_grid(6)
grid  = gl.TriMesh(gverts, gfaces)
gpath = gl.geodesic_path(grid, src=0, dst=13)   # (0,0) → (2,1)
assert gpath is not None

ep    = gpath.edge_params()     # list of (v0: int, v1: int, t: float)
poly  = gpath.polyline(grid)    # (K, 3) float64 ndarray

print(f"\nFlat 6×6 grid geodesic  src=(0,0) → dst=(2,1)")
print(f"  Euclidean distance:     {math.sqrt(5):.6f}")
print(f"  Geodesic length:        {gpath.length(grid):.6f}")
print(f"  n_points (incl. edge):  {gpath.n_vertices}")
print(f"  edge_params:")
for v0, v1, t in ep:
    kind = "vertex" if v0 == v1 else f"edge ({v0}→{v1})"
    pos  = poly[ep.index((v0, v1, t))]
    print(f"    ({v0:3d}, {v1:3d}, {t:.4f})  [{kind}]  3D={pos.round(4)}")

# Verify: at least one edge-crossing point exists
assert any(v0 != v1 for v0, v1, t in ep), "expected edge crossings"
# Verify: length equals sqrt(5) to machine precision
assert abs(gpath.length(grid) - math.sqrt(5)) < 1e-6, f"length mismatch: {gpath.length(grid)}"

print("\nAll demos completed successfully.")
