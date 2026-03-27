/// Geodesic path straightening via fan unfolding.
///
/// Paths are represented as sequences of `PathPoint`s — either mesh vertices
/// or points on mesh edge interiors (given as an edge + parameter t ∈ [0,1]).
///
/// The straightening algorithm (Polthier & Schmies 1998):
///   For each interior Vertex p_i with Vertex neighbours p_{i-1} and p_{i+1}:
///     1. Build the CCW (then CW) fan of triangles around p_i bounded by the
///        incoming and outgoing path edges.
///     2. Unfold the fan into the 2D plane.
///     3. Find where the straight line from p_{i-1} to p_{i+1} crosses a fan edge.
///     4. Replace p_i with that edge crossing (PathPoint::Edge).
///
/// This converges to a locally-shortest geodesic whose polyline passes through
/// edge interiors, not just along mesh edges.

use crate::mesh::{HalfedgeMesh, INVALID, dist3, sub3, dot3, norm3};

// ────────────────────────────────────────────────────────────────
// PathPoint
// ────────────────────────────────────────────────────────────────

/// A point on a geodesic: at a mesh vertex or on a mesh edge interior.
#[derive(Debug, Clone)]
pub enum PathPoint {
    /// At mesh vertex `v`.
    Vertex(usize),
    /// On the directed edge v0 → v1 at parameter t ∈ [0, 1]:
    ///   position = positions[v0] + t * (positions[v1] − positions[v0])
    Edge { v0: usize, v1: usize, t: f64 },
}

impl PathPoint {
    /// 3-D world position.
    pub fn position(&self, mesh: &HalfedgeMesh) -> [f64; 3] {
        match self {
            PathPoint::Vertex(v) => mesh.position(*v),
            PathPoint::Edge { v0, v1, t } => {
                let p0 = mesh.position(*v0);
                let p1 = mesh.position(*v1);
                [
                    p0[0] + t * (p1[0] - p0[0]),
                    p0[1] + t * (p1[1] - p0[1]),
                    p0[2] + t * (p1[2] - p0[2]),
                ]
            }
        }
    }

    /// Return the vertex index if this is a Vertex point.
    pub fn as_vertex(&self) -> Option<usize> {
        if let PathPoint::Vertex(v) = self { Some(*v) } else { None }
    }
}

// ────────────────────────────────────────────────────────────────
// GeodesicPath
// ────────────────────────────────────────────────────────────────

/// A geodesic path as a sequence of `PathPoint`s.
///
/// For a **closed** loop the last point is NOT repeated; the segment from
/// `points.last()` back to `points[0]` closes the loop implicitly.
#[derive(Debug, Clone)]
pub struct GeodesicPath {
    pub points: Vec<PathPoint>,
    pub is_closed: bool,
}

impl GeodesicPath {
    /// Construct an open path from a vertex sequence (e.g. from Dijkstra).
    pub fn open(vertices: Vec<usize>) -> Self {
        Self {
            points: vertices.into_iter().map(PathPoint::Vertex).collect(),
            is_closed: false,
        }
    }

    /// Construct a closed loop from a vertex sequence.
    pub fn closed(vertices: Vec<usize>) -> Self {
        Self {
            points: vertices.into_iter().map(PathPoint::Vertex).collect(),
            is_closed: true,
        }
    }

    pub fn len(&self) -> usize { self.points.len() }
    pub fn is_empty(&self) -> bool { self.points.is_empty() }

    /// Vertex indices only (Vertex points; Edge points are skipped).
    pub fn vertex_indices(&self) -> Vec<usize> {
        self.points.iter().filter_map(|p| p.as_vertex()).collect()
    }

    /// `(v0, v1, t)` for every point.  `Vertex(v)` maps to `(v, v, 0.0)`.
    pub fn edge_params(&self) -> Vec<(usize, usize, f64)> {
        self.points.iter().map(|p| match p {
            PathPoint::Vertex(v) => (*v, *v, 0.0),
            PathPoint::Edge { v0, v1, t } => (*v0, *v1, *t),
        }).collect()
    }

    /// 3-D polyline coordinates.  For closed paths the first point is
    /// appended at the end to close the polygon.
    pub fn to_polyline(&self, mesh: &HalfedgeMesh) -> Vec<[f64; 3]> {
        let mut pts: Vec<[f64; 3]> = self.points.iter()
            .map(|p| p.position(mesh)).collect();
        if self.is_closed && !pts.is_empty() {
            pts.push(pts[0]);
        }
        pts
    }

    /// Euclidean 3-D path length (sum of straight-line segment distances).
    pub fn euclidean_length(&self, mesh: &HalfedgeMesh) -> f64 {
        let pts: Vec<[f64; 3]> = self.points.iter()
            .map(|p| p.position(mesh)).collect();
        let n = pts.len();
        if n < 2 { return 0.0; }
        let mut total: f64 = (0..n - 1).map(|i| dist3(&pts[i], &pts[i + 1])).sum();
        if self.is_closed { total += dist3(&pts[n - 1], &pts[0]); }
        total
    }

    /// Metric length.  For pure-vertex paths the mesh edge weights are used;
    /// for paths with edge crossings the Euclidean length is returned
    /// (exact for the Euclidean metric, approximate otherwise).
    pub fn metric_length(&self, mesh: &HalfedgeMesh) -> f64 {
        if self.points.iter().all(|p| p.as_vertex().is_some()) {
            let verts = self.vertex_indices();
            let base = crate::geodesic::path_metric_length(mesh, &verts);
            let close = if self.is_closed && verts.len() >= 2 {
                let l = *verts.last().unwrap();
                let f = verts[0];
                if l == f { 0.0 } else { metric_edge_len_between(mesh, l, f) }
            } else {
                0.0
            };
            base + close
        } else {
            self.euclidean_length(mesh)
        }
    }
}

// ────────────────────────────────────────────────────────────────
// Internal mesh helpers
// ────────────────────────────────────────────────────────────────

fn metric_edge_len_between(mesh: &HalfedgeMesh, u: usize, v: usize) -> f64 {
    for h in mesh.outgoing_halfedges(u) {
        if mesh.dest(h) == v { return mesh.metric_edge_len(h); }
    }
    f64::INFINITY
}

/// Find the halfedge u → v, or INVALID.
fn find_halfedge(mesh: &HalfedgeMesh, u: usize, v: usize) -> usize {
    for h in mesh.outgoing_halfedges(u) {
        if mesh.dest(h) == v { return h; }
    }
    INVALID
}

// ────────────────────────────────────────────────────────────────
// 2-D geometry helpers
// ────────────────────────────────────────────────────────────────

#[inline] fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 { a[0]*b[1] - a[1]*b[0] }
#[inline] fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] { [a[0]-b[0], a[1]-b[1]] }

/// Intersect segment [p,q] with segment [s,t].
/// Returns `(r, u)` where r is the parameter along p→q and u along s→t.
fn seg_seg_intersect(
    p: [f64; 2], q: [f64; 2],
    s: [f64; 2], t: [f64; 2],
) -> Option<(f64, f64)> {
    let d = sub2(q, p);
    let e = sub2(t, s);
    let denom = cross2(d, e);
    if denom.abs() < 1e-14 { return None; }
    let f = sub2(s, p);
    Some((cross2(f, e) / denom, cross2(f, d) / denom))
}

// ────────────────────────────────────────────────────────────────
// Fan unfolding
// ────────────────────────────────────────────────────────────────

/// Unfold the **CCW** fan of triangles around vertex B = dest(h_in) = origin(h_out).
///
/// The fan sweeps CCW from `h_out` (B→C) to `twin(h_in)` (B→A).
///
/// Returns `(fan_he, verts_2d)` where:
/// * `fan_he[k]`   = k-th halfedge **from B** in fan order
/// * `verts_2d[k]` = 2-D position of `dest(fan_he[k])`
/// * `verts_2d[0]` = C placed at `(edge_len(h_out), 0)`
/// * `verts_2d.last()` = A placed at accumulated angle
fn unfold_fan_ccw(
    mesh: &HalfedgeMesh,
    h_in: usize,   // halfedge prev_v → B
    h_out: usize,  // halfedge B → next_v
) -> Option<(Vec<usize>, Vec<[f64; 2]>)> {
    let b = mesh.dest(h_in);
    if mesh.origin(h_out) != b { return None; }
    let h_target = mesh.twin(h_in); // B → A
    if h_target == INVALID { return None; }
    if mesh.is_boundary_he(h_out) { return None; }

    let mut fan_he: Vec<usize> = vec![h_out];
    let mut verts_2d: Vec<[f64; 2]> = vec![[mesh.edge_len(h_out), 0.0]];
    let mut cum = 0.0f64;
    let mut h = h_out;

    for _ in 0..512 {
        if mesh.is_boundary_he(h) { return None; }

        // Accumulate corner angle at B in face(h).
        // corner_angle(h) = angle at origin(h) = B ✓
        cum += mesh.corner_angle(h);

        // CCW step around B: twin(prev(h))
        let ph = mesh.prev(h);
        if ph == INVALID { return None; }
        let hn = mesh.twin(ph);
        if hn == INVALID { return None; }

        let r = mesh.edge_len(hn);
        verts_2d.push([r * cum.cos(), r * cum.sin()]);
        fan_he.push(hn);

        if hn == h_target { break; }
        h = hn;
        if fan_he.len() > 512 { return None; }
    }

    if fan_he.last().copied() != Some(h_target) || fan_he.len() < 2 {
        return None;
    }
    Some((fan_he, verts_2d))
}

/// Unfold the **CW** fan of triangles around B.
///
/// Same as CCW but traverses the opposite side; angles accumulate negatively
/// (vertices appear below the x-axis in the unfolded plane).
fn unfold_fan_cw(
    mesh: &HalfedgeMesh,
    h_in: usize,
    h_out: usize,
) -> Option<(Vec<usize>, Vec<[f64; 2]>)> {
    let b = mesh.dest(h_in);
    if mesh.origin(h_out) != b { return None; }
    let h_target = mesh.twin(h_in);
    if h_target == INVALID { return None; }

    let mut fan_he: Vec<usize> = vec![h_out];
    let mut verts_2d: Vec<[f64; 2]> = vec![[mesh.edge_len(h_out), 0.0]];
    let mut cum = 0.0f64;
    let mut h = h_out;

    for _ in 0..512 {
        // CW step around B: next(twin(h))
        let th = mesh.twin(h);
        if th == INVALID { return None; }
        if mesh.is_boundary_he(th) { return None; }
        let hn = mesh.next(th);
        if hn == INVALID { return None; }
        if mesh.origin(hn) != b { return None; }

        // corner_angle(hn) = angle at B in face(twin(h)) between B→dest(hn) and B→dest(h)
        cum -= mesh.corner_angle(hn);

        let r = mesh.edge_len(hn);
        verts_2d.push([r * cum.cos(), r * cum.sin()]);
        fan_he.push(hn);

        if hn == h_target { break; }
        h = hn;
        if fan_he.len() > 512 { return None; }
    }

    if fan_he.last().copied() != Some(h_target) || fan_he.len() < 2 {
        return None;
    }
    Some((fan_he, verts_2d))
}

/// Given an unfolded fan, find all crossings of the straight line from
/// `verts_2d.last()` (A) to `verts_2d[0]` (C) with the **radial** fan edges.
///
/// Radial edges: from B=(0,0) to `verts_2d[k]` for k in `1..n-2`.
/// (k=0 is B→C and k=n-1 is B→A — these are the path edges, not interior.)
///
/// Each crossing at segment-parameter `t` along B→V_k becomes a
/// `PathPoint::Edge { v0: B, v1: V_k, t }`.  Multiple crossings are
/// returned sorted by `ray_s` (order A→C along the geodesic).
fn find_fan_crossings(
    mesh: &HalfedgeMesh,
    fan_he: &[usize],
    verts_2d: &[[f64; 2]],
) -> Vec<PathPoint> {
    let n = verts_2d.len();
    if n < 3 { return Vec::new(); }

    let a2 = verts_2d[n - 1]; // A in 2D
    let c2 = verts_2d[0];     // C in 2D
    let b2 = [0.0f64, 0.0];   // B at origin

    const EPS: f64 = 1e-9;
    let mut hits: Vec<(f64, PathPoint)> = Vec::new();

    // Check each radial edge B → verts_2d[k] for interior k
    for k in 1..n - 1 {
        let vk = verts_2d[k];

        if let Some((ray_s, seg_t)) = seg_seg_intersect(a2, c2, b2, vk) {
            if ray_s > EPS && ray_s < 1.0 - EPS && seg_t > EPS && seg_t < 1.0 - EPS {
                // fan_he[k] is the halfedge B → V_k
                let h = fan_he[k];
                let v0 = mesh.origin(h); // B
                let v1 = mesh.dest(h);   // V_k
                hits.push((ray_s, PathPoint::Edge { v0, v1, t: seg_t }));
            }
        }
    }

    // Sort by ray_s so crossings are ordered A → C
    hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    hits.into_iter().map(|(_, pt)| pt).collect()
}

// ────────────────────────────────────────────────────────────────
// Per-vertex straightening
// ────────────────────────────────────────────────────────────────

/// Try to straighten the path at Vertex `v` given its neighbours.
///
/// Returns the replacement PathPoints (0 = already straight, 1+ = edge crossings).
/// Currently only handles Vertex neighbours (typical after Dijkstra).
fn straighten_at_vertex(
    mesh: &HalfedgeMesh,
    prev_pt: &PathPoint,
    v: usize,
    next_pt: &PathPoint,
) -> Vec<PathPoint> {
    let prev_v = match prev_pt.as_vertex() { Some(u) => u, None => return Vec::new() };
    let next_v = match next_pt.as_vertex() { Some(u) => u, None => return Vec::new() };
    if prev_v == v || next_v == v || prev_v == next_v { return Vec::new(); }

    let h_in  = find_halfedge(mesh, prev_v, v);
    let h_out = find_halfedge(mesh, v, next_v);
    if h_in == INVALID || h_out == INVALID { return Vec::new(); }

    // Try CCW fan first (left turn at v)
    if let Some((fan_he, verts_2d)) = unfold_fan_ccw(mesh, h_in, h_out) {
        let pts = find_fan_crossings(mesh, &fan_he, &verts_2d);
        if !pts.is_empty() { return pts; }
    }

    // Try CW fan (right turn at v)
    if let Some((fan_he, verts_2d)) = unfold_fan_cw(mesh, h_in, h_out) {
        let pts = find_fan_crossings(mesh, &fan_he, &verts_2d);
        if !pts.is_empty() { return pts; }
    }

    Vec::new() // v is already locally straight (or a cone point)
}

// ────────────────────────────────────────────────────────────────
// Configuration
// ────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FlipConfig {
    pub max_iterations: usize,
    pub rel_tol: f64,
    /// Kept for API compatibility; unused in the unfolding algorithm.
    pub angle_eps: f64,
}

impl Default for FlipConfig {
    fn default() -> Self {
        Self { max_iterations: 10_000, rel_tol: 1e-9, angle_eps: 1e-6 }
    }
}

// ────────────────────────────────────────────────────────────────
// Main shortening loop
// ────────────────────────────────────────────────────────────────

/// Shorten a geodesic path by iteratively replacing interior Vertex points
/// with edge-crossing PathPoints via fan unfolding.
pub fn shorten_path(
    mesh: &HalfedgeMesh,
    path: &GeodesicPath,
    config: &FlipConfig,
) -> GeodesicPath {
    let mut pts = path.points.clone();
    let is_closed = path.is_closed;

    if pts.len() < 3 {
        return path.clone();
    }

    let mut prev_len = f64::INFINITY;

    for _iter in 0..config.max_iterations {
        let n = pts.len();
        if n < 3 { break; }

        let mut changed = false;

        let mut i = 0;
        while i < pts.len() {
            let n_cur = pts.len();
            if n_cur < 3 { break; }

            // Determine indices, accounting for closed paths
            let (prev_i, cur_i, next_i) = if is_closed {
                let c = i % n_cur;
                let p = if c == 0 { n_cur - 1 } else { c - 1 };
                let nx = if c == n_cur - 1 { 0 } else { c + 1 };
                (p, c, nx)
            } else {
                if i == 0 || i >= n_cur - 1 { i += 1; continue; }
                (i - 1, i, i + 1)
            };

            // Only straighten at Vertex points
            let v = match &pts[cur_i] {
                PathPoint::Vertex(v) => *v,
                _ => { i += 1; continue; }
            };

            let prev_pt = pts[prev_i].clone();
            let next_pt = pts[next_i].clone();

            let replacements = straighten_at_vertex(mesh, &prev_pt, v, &next_pt);
            if !replacements.is_empty() {
                let rep_len = replacements.len();
                pts.splice(cur_i..=cur_i, replacements.into_iter());
                changed = true;
                // Skip the newly-inserted edge points (they aren't Vertex points)
                i = cur_i + rep_len;
                continue;
            }

            i += 1;
        }

        let cur_len = GeodesicPath { points: pts.clone(), is_closed }
            .euclidean_length(mesh);

        if cur_len < 1e-15 { break; }
        let rel_dec = (prev_len - cur_len) / prev_len.max(cur_len);
        if !changed || rel_dec.abs() < config.rel_tol { break; }
        prev_len = cur_len;
    }

    GeodesicPath { points: pts, is_closed }
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

// ────────────────────────────────────────────────────────────────
// High-level API
// ────────────────────────────────────────────────────────────────

pub fn geodesic_path(
    mesh: &HalfedgeMesh,
    src: usize,
    dst: usize,
    config: &FlipConfig,
) -> Option<GeodesicPath> {
    let vpath = crate::geodesic::shortest_path(mesh, src, dst)?;
    Some(shorten_path(mesh, &GeodesicPath::open(vpath), config))
}

pub fn geodesic_loop(
    mesh: &HalfedgeMesh,
    seed: usize,
    config: &FlipConfig,
) -> Option<GeodesicPath> {
    let vloop = crate::geodesic::shortest_loop_through(mesh, seed)?;
    Some(shorten_loop(mesh, &GeodesicPath::closed(vloop), config))
}

pub fn geodesic_loop_with_cut(
    mesh: &HalfedgeMesh,
    seed: usize,
    cut_edges: &std::collections::HashSet<usize>,
    config: &FlipConfig,
) -> Option<GeodesicPath> {
    let vloop = crate::geodesic::shortest_loop_crossing_cut(mesh, seed, cut_edges)?;
    Some(shorten_loop(mesh, &GeodesicPath::closed(vloop), config))
}

// ────────────────────────────────────────────────────────────────
// Statistics (kept for Python API compatibility)
// ────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct GeodesicStats {
    pub length: f64,
    pub n_points: usize,
    pub is_closed: bool,
    /// Min angle at vertex points (radians). f64::INFINITY if no vertex turns.
    pub min_angle: f64,
    /// Max angle at vertex points (radians).
    pub max_angle: f64,
}

pub fn compute_stats(mesh: &HalfedgeMesh, path: &GeodesicPath) -> GeodesicStats {
    let length = path.metric_length(mesh);
    let n = path.points.len();
    let mut min_angle = f64::INFINITY;
    let mut max_angle = 0.0f64;

    let range: Box<dyn Iterator<Item = usize>> = if path.is_closed {
        Box::new(0..n)
    } else {
        Box::new(1..n.saturating_sub(1))
    };

    for i in range {
        let prev_i = if i == 0 { n - 1 } else { i - 1 };
        let next_i = if i == n - 1 { 0 } else { i + 1 };
        if let (Some(pv), Some(v), Some(nv)) = (
            path.points[prev_i].as_vertex(),
            path.points[i].as_vertex(),
            path.points[next_i].as_vertex(),
        ) {
            let (l, r) = wedge_angles(mesh, pv, v, nv);
            let ma = l.min(r);
            min_angle = min_angle.min(ma);
            max_angle = max_angle.max(ma);
        }
    }

    GeodesicStats { length, n_points: n, is_closed: path.is_closed, min_angle, max_angle }
}

// ────────────────────────────────────────────────────────────────
// Angle utilities (kept for stats / backward compat)
// ────────────────────────────────────────────────────────────────

pub fn wedge_angles(
    mesh: &HalfedgeMesh,
    prev_v: usize, v: usize, next_v: usize,
) -> (f64, f64) {
    let pv = mesh.position(v);
    let pp = mesh.position(prev_v);
    let pn = mesh.position(next_v);
    let dp = sub3(&pp, &pv);
    let dn = sub3(&pn, &pv);
    let cos_a = dot3(&dp, &dn) / (norm3(&dp) * norm3(&dn) + 1e-30);
    let angle_between = cos_a.clamp(-1.0, 1.0).acos();
    let angle_sum = mesh.vertex_angle_sum(v);
    (angle_between, (angle_sum - angle_between).abs())
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TurnType { Straight, LeftTurn, RightTurn }

pub fn classify_turn(
    mesh: &HalfedgeMesh,
    prev_v: usize, v: usize, next_v: usize,
    eps: f64,
) -> (TurnType, f64) {
    let (left, right) = wedge_angles(mesh, prev_v, v, next_v);
    let min_angle = left.min(right);
    if min_angle >= std::f64::consts::PI - eps { (TurnType::Straight, min_angle) }
    else if left < right { (TurnType::LeftTurn, left) }
    else { (TurnType::RightTurn, right) }
}

pub fn is_delaunay_edge(mesh: &HalfedgeMesh, h: usize) -> bool {
    if mesh.is_boundary_edge(h) { return true; }
    let t = mesh.twin(h);
    let alpha = mesh.corner_angle(mesh.next(h));
    let beta  = mesh.corner_angle(mesh.next(t));
    alpha + beta <= std::f64::consts::PI + 1e-10
}
