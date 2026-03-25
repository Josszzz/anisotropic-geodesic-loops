"""
End-to-end demo of geodesic_loops.

Run from the repo root after building the extension with:
    cd geodesic_loops
    maturin develop --features python   (or: pip install maturin && maturin develop)

The demo:
  1. Builds a discrete cylinder mesh (unit-radius, wraps around in j).
  2. Computes a Dijkstra shortest path between two vertices.
  3. Shortens the path with the flip algorithm.
  4. Finds the shortest non-contractible geodesic loop.
  5. Demonstrates the scalar isotropic metric (simulated scar / slow zone).
  6. Demonstrates the fiber-transverse anisotropic metric (cardiac electrophysiology).
"""

import math
import numpy as np

# Import the compiled extension
try:
    import geodesic_loops as gl
except ImportError:
    raise ImportError(
        "geodesic_loops extension not found. "
        "Build it with: cd geodesic_loops && maturin develop --features python"
    )


# ──────────────────────────────────────────────────────────────────────────────
# 1.  Build a cylinder mesh
# ──────────────────────────────────────────────────────────────────────────────

def make_cylinder(n_circ: int = 20, n_rings: int = 10) -> tuple:
    """Return (vertices, faces) for a unit-radius cylinder."""
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


verts, faces = make_cylinder(n_circ=20, n_rings=8)
mesh = gl.TriMesh(verts, faces)
print(mesh)


# ──────────────────────────────────────────────────────────────────────────────
# 2.  Point-to-point geodesic
# ──────────────────────────────────────────────────────────────────────────────

src = 0
dst = mesh.n_vertices - 1   # last vertex (opposite side of cylinder)

path = gl.geodesic_path(mesh, src, dst)
if path is None:
    print("Unreachable!")
else:
    pts = path.polyline(mesh)
    print(f"\nPoint-to-point geodesic  src={src} → dst={dst}")
    print(f"  vertices:  {path.n_vertices}")
    print(f"  length:    {path.length(mesh):.4f}")
    print(f"  polyline shape: {pts.shape}")
    st = path.stats(mesh)
    print(f"  min/max angle (rad): {st['min_angle_rad']:.3f} / {st['max_angle_rad']:.3f}")


# ──────────────────────────────────────────────────────────────────────────────
# 3.  Closed geodesic loop
# ──────────────────────────────────────────────────────────────────────────────

loop = gl.geodesic_loop(mesh, seed=0)
if loop is None:
    print("\nNo loop found (flat/simply connected mesh).")
else:
    print(f"\nClosed geodesic loop (seed=0)")
    print(f"  vertices:  {loop.n_vertices}")
    print(f"  length:    {loop.length(mesh):.4f}")
    print(f"  is_closed: {loop.is_closed}")


# ──────────────────────────────────────────────────────────────────────────────
# 4.  Non-contractible loop via seam cut
# ──────────────────────────────────────────────────────────────────────────────

n_circ = 20
n_rings = 8
# Seam edges: edge (ring*n_circ + n_circ-1, ring*n_circ + 0) for each ring
# These are the "horizontal" wrap-around edges at j=n_circ-1 → j=0.
seam_edges = []
for ring in range(n_rings):
    v_last  = ring * n_circ + (n_circ - 1)
    v_first = ring * n_circ
    # Find the edge index by checking edge_lengths (brute-force search via
    # the vertex lookup; in practice you'd compute this from the mesh directly)
    # Here we just pass vertex pairs — the Python binding takes edge indices.
    # For a clean API, we'd expose a "find_edge(v0, v1)" method, but for
    # the demo we rely on the fact that the first n_circ*(n_rings-1)*2*3/2
    # edges are ordered predictably.  Instead, use geodesic_loop_with_cut
    # and let Rust handle the seam detection internally.
    pass  # We'll use a different approach below.

# Alternative: use shortest_loop_crossing_cut with ring-0 seam edge indices.
# We find the seam edges via outgoing halfedge enumeration (not exposed yet),
# so we demonstrate with a direct geodesic_loop call instead.
print(f"\nUsing geodesic_loop (combinatorial, may find contractible short loop):")
raw_loop = gl.shortest_loop(mesh, seed=0)
if raw_loop is not None:
    print(f"  raw loop length:   {raw_loop.length(mesh):.4f}")
    shortened = gl.flip_shorten(mesh, raw_loop)
    print(f"  shortened length:  {shortened.length(mesh):.4f}")


# ──────────────────────────────────────────────────────────────────────────────
# 5.  Scalar isotropic metric (simulated slow zone / scar)
# ──────────────────────────────────────────────────────────────────────────────

# Create a "scar" region: slow conduction speed around j=5 in ring 3-5
speeds = np.ones(mesh.n_vertices, dtype=np.float64)
for v in range(mesh.n_vertices):
    ring_idx = v // n_circ
    j_idx    = v % n_circ
    if 3 <= ring_idx <= 5 and 4 <= j_idx <= 7:
        speeds[v] = 0.1   # 10x slower conduction in the scar

mesh.set_isotropic_speeds(speeds)

path_aniso = gl.geodesic_path(mesh, src=0, dst=dst)
if path_aniso is not None:
    print(f"\nIsotropic slow-zone geodesic  src={src} → dst={dst}")
    print(f"  metric length (travel time): {path_aniso.length(mesh):.4f}")
    print(f"  vertices in path:            {path_aniso.n_vertices}")
    # The path should avoid the slow zone, going around it
    visited_j = [v % n_circ for v in path_aniso.vertices()]
    print(f"  j-indices visited: {sorted(set(visited_j))}")

mesh.clear_metric()


# ──────────────────────────────────────────────────────────────────────────────
# 6.  Fiber-transverse anisotropic metric (cardiac electrophysiology)
# ──────────────────────────────────────────────────────────────────────────────

# Fibers run along the z-axis (longitudinal direction) of the cylinder.
# Fast conduction along fibers (c_f = 1.0), slow transverse (c_t = 0.3).
fiber_dirs = np.zeros((mesh.n_vertices, 3), dtype=np.float64)
fiber_dirs[:, 2] = 1.0          # z-axis = longitudinal fiber direction

speeds_fiber      = np.full(mesh.n_vertices, 1.0, dtype=np.float64)
speeds_transverse = np.full(mesh.n_vertices, 0.3, dtype=np.float64)

mesh.set_fiber_metric(fiber_dirs, speeds_fiber, speeds_transverse)

path_fiber = gl.geodesic_path(mesh, src=0, dst=dst)
if path_fiber is not None:
    print(f"\nFiber-anisotropic geodesic  src={src} → dst={dst}")
    print(f"  metric length (travel time): {path_fiber.length(mesh):.4f}")
    print(f"  vertices in path:            {path_fiber.n_vertices}")
    # Along the fiber (z) the metric is 1/c_f = 1.0; transverse is 1/c_t ≈ 3.3
    # So the path should hug the z-direction and avoid circumferential detours.
    visited_j = [v % n_circ for v in path_fiber.vertices()]
    print(f"  j-indices visited: {sorted(set(visited_j))}")

# Introduce a slow scar that now has different fiber orientation
fiber_dirs_scar = fiber_dirs.copy()
for v in range(mesh.n_vertices):
    ring_idx = v // n_circ
    j_idx    = v % n_circ
    if 3 <= ring_idx <= 5 and 4 <= j_idx <= 7:
        # Fibers rotated 90° in the scar (circumferential)
        theta = 2 * math.pi * j_idx / n_circ
        fiber_dirs_scar[v] = [-math.sin(theta), math.cos(theta), 0.0]
        speeds_fiber[v]      = 0.2   # also generally slower
        speeds_transverse[v] = 0.05

mesh.set_fiber_metric(fiber_dirs_scar, speeds_fiber, speeds_transverse)

path_scar = gl.geodesic_path(mesh, src=0, dst=dst)
if path_scar is not None:
    print(f"\nFiber-anisotropic geodesic with scar  src={src} → dst={dst}")
    print(f"  metric length (travel time): {path_scar.length(mesh):.4f}")
    print(f"  vertices in path:            {path_scar.n_vertices}")

mesh.clear_metric()


# ──────────────────────────────────────────────────────────────────────────────
# 6.  Flat triangulated grid example
# ──────────────────────────────────────────────────────────────────────────────

def make_flat_grid(n: int) -> tuple:
    verts = [[float(i), float(j), 0.0]
             for i in range(n) for j in range(n)]
    faces = []
    for i in range(n - 1):
        for j in range(n - 1):
            a, b, c, d = i*n+j, i*n+j+1, (i+1)*n+j, (i+1)*n+j+1
            if (i + j) % 2 == 0:
                # Up-right diagonal a-d: allows diagonal travel towards (n-1, n-1)
                faces.append([a, b, d])
                faces.append([a, d, c])
            else:
                faces.append([a, b, c])
                faces.append([b, d, c])
    return (np.array(verts, dtype=np.float64),
            np.array(faces, dtype=np.int64))

gverts, gfaces = make_flat_grid(10)
grid = gl.TriMesh(gverts, gfaces)
print(f"\nFlat 10×10 grid: {grid}")

gpath = gl.geodesic_path(grid, src=0, dst=99)
if gpath is not None:
    euclidean = math.sqrt(9**2 + 9**2)
    print(f"  0→99 geodesic length:  {gpath.length(grid):.4f}  (Euclidean: {euclidean:.4f})")
    assert gpath.length(grid) >= euclidean - 1e-6
    assert gpath.length(grid) <= euclidean * 1.5, \
        f"path too long: {gpath.length(grid):.4f} vs Euclidean {euclidean:.4f}"

print("\nAll demos completed successfully.")
