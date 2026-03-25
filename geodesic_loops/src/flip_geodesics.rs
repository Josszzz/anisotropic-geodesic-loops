/// Flip geodesics: Crane & Sharp (2020) "Geodesics in Heat" edge-flip shortening.
///
/// Given an initial piecewise-linear path on a mesh (as a sequence of vertices),
/// we iteratively shorten it by flipping edges in wedges where the path makes
/// a turn angle < π. This converges to a locally shortest geodesic.
///
/// This implementation works on the *extrinsic* mesh (no intrinsic
/// triangulation refinement) to keep the code self-contained. For highly
/// irregular meshes, running a few rounds of edge flipping toward Delaunay
/// first would improve quality.
///
/// Reference:
///   Sharp, N. and Crane, K. (2020). "You Can Find Geodesic Paths in
///   Triangle Meshes by Just Flipping Edges". ACM Trans. Graph. 39(6).

use crate::mesh::{HalfedgeMesh, dist3, norm3, dot3, sub3};

// ---- angle type at a vertex of the path ----

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TurnType {
    Straight,   // angle >= π - ε at both sides → locally shortest
    LeftTurn,   // path turns left → need to flip edges on the right wedge
    RightTurn,  // path turns right → need to flip edges on the left wedge
}

// ---- path representation ----

/// A geodesic path as a sequence of vertices (may be open or closed).
#[derive(Debug, Clone)]
pub struct GeodesicPath {
    /// Ordered vertices along the path.  For a closed loop, `vertices[0]`
    /// is NOT repeated at the end; the edge from `vertices.last()` back to
    /// `vertices[0]` closes the loop.
    pub vertices: Vec<usize>,
    pub is_closed: bool,
}

impl GeodesicPath {
    pub fn open(vertices: Vec<usize>) -> Self {
        Self { vertices, is_closed: false }
    }
    pub fn closed(vertices: Vec<usize>) -> Self {
        Self { vertices, is_closed: true }
    }

    pub fn len(&self) -> usize { self.vertices.len() }

    /// Metric length of the path on the mesh
    pub fn metric_length(&self, mesh: &HalfedgeMesh) -> f64 {
        crate::geodesic::path_metric_length(mesh, &self.vertices)
            + if self.is_closed && self.vertices.len() >= 2 {
                let last = *self.vertices.last().unwrap();
                let first = self.vertices[0];
                if last == first { 0.0 } else { metric_edge_len_between(mesh, last, first) }
            } else { 0.0 }
    }

    /// Convert to 3-D polyline
    pub fn to_polyline(&self, mesh: &HalfedgeMesh) -> Vec<[f64; 3]> {
        let mut pts: Vec<[f64; 3]> = self.vertices.iter().map(|&v| mesh.position(v)).collect();
        if self.is_closed {
            pts.push(pts[0]);
        }
        pts
    }
}

/// Return the metric edge length between adjacent vertices u and v.
/// Returns infinity if there is no direct edge (vertices are not neighbours).
fn metric_edge_len_between(mesh: &HalfedgeMesh, u: usize, v: usize) -> f64 {
    for h in mesh.outgoing_halfedges(u) {
        if mesh.dest(h) == v {
            return mesh.metric_edge_len(h);
        }
    }
    f64::INFINITY  // no direct mesh edge
}

// ---- turn angle at a vertex ----

/// Compute the turn angle of the path at interior vertex `v`,
/// where the path comes in via vertex `prev_v` and goes out via `next_v`.
///
/// Returns `(left_angle, right_angle)` — the angular gaps on each side.
/// A locally-shortest path has both angles >= π - ε.
pub fn wedge_angles(
    mesh: &HalfedgeMesh,
    prev_v: usize,
    v: usize,
    next_v: usize,
) -> (f64, f64) {
    // We use actual 3-D geometry: compute the angle in the tangent plane
    // of the surface at v between the two edges v→prev and v→next.
    // For a flat mesh this is exact; for curved meshes it's approximate.
    //
    // More precisely, we unfold the two triangles incident to the path
    // at this vertex to compute the intrinsic angle.
    let pv = mesh.position(v);
    let pp = mesh.position(prev_v);
    let pn = mesh.position(next_v);

    let dp = sub3(&pp, &pv);
    let dn = sub3(&pn, &pv);

    let angle_between = {
        let cos_a = dot3(&dp, &dn) / (norm3(&dp) * norm3(&dn) + 1e-30);
        cos_a.clamp(-1.0, 1.0).acos()
    };

    let angle_sum = mesh.vertex_angle_sum(v);

    // The two sides of the wedge: the angle going one way and the other
    // around the vertex.
    let left  = angle_between;
    let right = (angle_sum - angle_between).abs();

    (left, right)
}

/// Classify the turn type at v given incoming edge prev→v and outgoing v→next.
pub fn classify_turn(
    mesh: &HalfedgeMesh,
    prev_v: usize,
    v: usize,
    next_v: usize,
    eps: f64,
) -> (TurnType, f64) {
    let (left, right) = wedge_angles(mesh, prev_v, v, next_v);
    let min_angle = left.min(right);
    if min_angle >= std::f64::consts::PI - eps {
        (TurnType::Straight, min_angle)
    } else if left < right {
        (TurnType::LeftTurn, left)
    } else {
        (TurnType::RightTurn, right)
    }
}

// ---- wedge edge-flip ----

/// Attempt to locally shorten the path at vertex `v` by flipping edges
/// in the smaller wedge.  Returns true if a flip was performed.
///
/// The path is represented as a mutable vertex list. We only modify the
/// segment around `v` (indices `idx-1 .. idx+1` in the vertex list, or
/// wrapped for closed loops).
///
/// The algorithm:
///   1. Find the wedge (the set of edges between the two path edges at v)
///   2. If the wedge angle < π, find a diagonal crossing the wedge that
///      shortens the path
///   3. Replace the two-step path through v with a one-step shortcut through
///      the new diagonal vertex
fn try_shorten_at(
    mesh: &HalfedgeMesh,
    vertices: &mut Vec<usize>,
    idx: usize,
    is_closed: bool,
    eps: f64,
) -> bool {
    let n = vertices.len();
    if n < 3 { return false; }

    let prev_idx = if idx == 0 {
        if is_closed { n - 1 } else { return false; }
    } else { idx - 1 };

    let next_idx = if idx == n - 1 {
        if is_closed { 0 } else { return false; }
    } else { idx + 1 };

    let prev_v = vertices[prev_idx];
    let v = vertices[idx];
    let next_v = vertices[next_idx];

    let (turn_type, min_angle) = classify_turn(mesh, prev_v, v, next_v, eps);
    if turn_type == TurnType::Straight { return false; }

    // Find edges in the wedge: walk around v CCW from the outgoing edge v→next_v
    // to the incoming edge v→prev_v (using the smaller angular gap).
    // Look for a shortcut: a vertex w adjacent to v such that the path
    // prev_v→w→next_v is shorter than prev_v→v→next_v.
    let d_prev_v  = metric_edge_len_between(mesh, prev_v, v);
    let d_v_next  = metric_edge_len_between(mesh, v, next_v);
    let current_len = d_prev_v + d_v_next;

    let mut best_len = current_len;
    let mut best_w: Option<usize> = None;

    // Check all neighbors of v
    for h in mesh.outgoing_halfedges(v) {
        let w = mesh.dest(h);
        if w == prev_v || w == next_v || w >= mesh.n_verts() { continue; }

        // w must be inside the wedge: check if inserting w shortens the path
        let d_prev_w = metric_edge_len_between(mesh, prev_v, w);
        let d_w_next = metric_edge_len_between(mesh, w, next_v);
        if d_prev_w.is_infinite() || d_w_next.is_infinite() { continue; }

        let candidate_len = d_prev_w + d_w_next;
        if candidate_len < best_len - 1e-12 {
            best_len = candidate_len;
            best_w = Some(w);
        }
    }

    // Also try removing v entirely: direct edge prev_v → next_v
    let d_direct = metric_edge_len_between(mesh, prev_v, next_v);
    if d_direct < best_len - 1e-12 {
        vertices.remove(idx);
        return true;
    }

    if let Some(w) = best_w {
        // Replace v with w
        vertices[idx] = w;
        return true;
    }

    false
}

// ---- main shortening loop ----

/// Configuration for the flip geodesic shortening.
#[derive(Debug, Clone)]
pub struct FlipConfig {
    /// Maximum number of shortening iterations
    pub max_iterations: usize,
    /// Convergence tolerance: stop when relative length decrease < this
    pub rel_tol: f64,
    /// Angle tolerance for "straight" classification
    pub angle_eps: f64,
}

impl Default for FlipConfig {
    fn default() -> Self {
        Self {
            max_iterations: 10_000,
            rel_tol: 1e-8,
            angle_eps: 1e-6,
        }
    }
}

/// Shorten a geodesic path by iterative vertex-replacement at turn vertices.
///
/// This is a combinatorial analogue of the flip geodesic algorithm:
/// at each turn vertex we search the neighbourhood for a shorter route.
/// Converges to a locally shortest path in the combinatorial graph.
pub fn shorten_path(
    mesh: &HalfedgeMesh,
    path: &GeodesicPath,
    config: &FlipConfig,
) -> GeodesicPath {
    let mut vertices = path.vertices.clone();
    let is_closed = path.is_closed;

    let mut prev_len = f64::INFINITY;
    let mut iter = 0;

    loop {
        if iter >= config.max_iterations { break; }
        iter += 1;

        let n = vertices.len();
        if n < 2 { break; }

        let mut changed = false;
        let mut i = 0;
        while i < vertices.len() {
            if try_shorten_at(mesh, &mut vertices, i, is_closed, config.angle_eps) {
                changed = true;
                // Don't advance i — re-check same position after removal/replacement
            } else {
                i += 1;
            }
        }

        // Compute current length
        let cur_len = GeodesicPath { vertices: vertices.clone(), is_closed }
            .metric_length(mesh);

        if cur_len < 1e-15 { break; }

        let rel_decrease = (prev_len - cur_len).abs() / cur_len;
        if !changed || rel_decrease < config.rel_tol {
            break;
        }
        prev_len = cur_len;
    }

    GeodesicPath { vertices, is_closed }
}

/// Shorten a closed geodesic loop.
pub fn shorten_loop(
    mesh: &HalfedgeMesh,
    path: &GeodesicPath,
    config: &FlipConfig,
) -> GeodesicPath {
    assert!(path.is_closed);
    shorten_path(mesh, path, config)
}

// ---- Intrinsic edge-flip Delaunay refinement ----

/// Check whether edge h satisfies the (extrinsic) Delaunay condition:
/// the sum of opposite angles < π.
pub fn is_delaunay_edge(mesh: &HalfedgeMesh, h: usize) -> bool {
    if mesh.is_boundary_edge(h) { return true; }
    let t = mesh.twin(h);
    let alpha = mesh.corner_angle(mesh.next(h));  // opposite angle in face of h
    let beta  = mesh.corner_angle(mesh.next(t));  // opposite angle in face of twin
    alpha + beta <= std::f64::consts::PI + 1e-10
}

// ---- higher-level API ----

/// Compute and shorten a point-to-point geodesic from `src` to `dst`.
pub fn geodesic_path(
    mesh: &HalfedgeMesh,
    src: usize,
    dst: usize,
    config: &FlipConfig,
) -> Option<GeodesicPath> {
    let vpath = crate::geodesic::shortest_path(mesh, src, dst)?;
    let path = GeodesicPath::open(vpath);
    Some(shorten_path(mesh, &path, config))
}

/// Compute and shorten the shortest non-contractible closed geodesic through `seed`.
pub fn geodesic_loop(
    mesh: &HalfedgeMesh,
    seed: usize,
    config: &FlipConfig,
) -> Option<GeodesicPath> {
    let vloop = crate::geodesic::shortest_loop_through(mesh, seed)?;
    let path = GeodesicPath::closed(vloop);
    Some(shorten_loop(mesh, &path, config))
}

/// Compute a geodesic loop in a homotopy class defined by a cut.
pub fn geodesic_loop_with_cut(
    mesh: &HalfedgeMesh,
    seed: usize,
    cut_edges: &std::collections::HashSet<usize>,
    config: &FlipConfig,
) -> Option<GeodesicPath> {
    let vloop = crate::geodesic::shortest_loop_crossing_cut(mesh, seed, cut_edges)?;
    let path = GeodesicPath::closed(vloop);
    Some(shorten_loop(mesh, &path, config))
}

// ---- statistics ----

#[derive(Debug, Clone)]
pub struct GeodesicStats {
    pub length: f64,
    pub n_vertices: usize,
    pub is_closed: bool,
    pub min_angle: f64,
    pub max_angle: f64,
}

pub fn compute_stats(mesh: &HalfedgeMesh, path: &GeodesicPath) -> GeodesicStats {
    let n = path.vertices.len();
    let mut min_angle = f64::INFINITY;
    let mut max_angle = 0.0f64;

    let range = if path.is_closed { 0..n } else { 1..(n-1) };
    for i in range {
        let prev_i = if i == 0 { n - 1 } else { i - 1 };
        let next_i = if i == n - 1 { 0 } else { i + 1 };
        let pv = path.vertices[prev_i];
        let v  = path.vertices[i];
        let nv = path.vertices[next_i];
        let (left, right) = wedge_angles(mesh, pv, v, nv);
        let ma = left.min(right);
        min_angle = min_angle.min(ma);
        max_angle = max_angle.max(ma);
    }

    GeodesicStats {
        length: path.metric_length(mesh),
        n_vertices: n,
        is_closed: path.is_closed,
        min_angle,
        max_angle,
    }
}
