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

    /// Metric length of the path.
    ///
    /// For each segment, the metric length is approximated as:
    ///   Euclidean_segment_length × (metric_edge_len / euclidean_edge_len)
    /// averaged over the mesh edges that the segment's endpoints lie on.
    /// For the Euclidean metric this reduces to the Euclidean length exactly.
    pub fn metric_length(&self, mesh: &HalfedgeMesh) -> f64 {
        // For pure-vertex paths, use the exact mesh edge weights.
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
            return base + close;
        }
        // For paths with edge crossings, sum metric-weighted segment lengths.
        let positions: Vec<[f64; 3]> = self.points.iter()
            .map(|p| p.position(mesh)).collect();
        let n = positions.len();
        if n < 2 { return 0.0; }
        let seg_count = if self.is_closed { n } else { n - 1 };
        let mut total = 0.0;
        for s in 0..seg_count {
            let j = (s + 1) % n;
            let euc_d = dist3(&positions[s], &positions[j]);
            if euc_d < 1e-30 { continue; }
            // Compute the metric/euclidean ratio at each endpoint
            let r0 = metric_ratio_at(&self.points[s], mesh);
            let r1 = metric_ratio_at(&self.points[j], mesh);
            // Average ratio (trapezoidal approximation)
            total += euc_d * 0.5 * (r0 + r1);
        }
        total
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

/// Compute the metric/euclidean ratio at a path point.
/// For Vertex(v), uses the average ratio of incident edges.
/// For Edge{v0,v1,t}, interpolates the ratio of that edge.
fn metric_ratio_at(pt: &PathPoint, mesh: &HalfedgeMesh) -> f64 {
    match pt {
        PathPoint::Vertex(v) => {
            let mut sum_m = 0.0;
            let mut sum_e = 0.0;
            for h in mesh.outgoing_halfedges(*v) {
                let e = mesh.edge_len(h);
                let m = mesh.metric_edge_len(h);
                sum_e += e;
                sum_m += m;
            }
            if sum_e > 1e-30 { sum_m / sum_e } else { 1.0 }
        }
        PathPoint::Edge { v0, v1, .. } => {
            for h in mesh.outgoing_halfedges(*v0) {
                if mesh.dest(h) == *v1 {
                    let e = mesh.edge_len(h);
                    return if e > 1e-30 { mesh.metric_edge_len(h) / e } else { 1.0 };
                }
            }
            1.0
        }
    }
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

/// Unfold the **CCW** fan of triangles around vertex B.
///
/// The fan sweeps CCW from `h_out` (B→C) to `h_target` (B→A).
///
/// Returns `(fan_he, verts_2d)` where:
/// * `fan_he[k]`   = k-th halfedge **from B** in fan order
/// * `verts_2d[k]` = 2-D position of `dest(fan_he[k])`
/// * `verts_2d[0]` = C placed at `(edge_len(h_out), 0)`
/// * `verts_2d.last()` = dest(h_target) placed at accumulated angle
fn unfold_fan_ccw(
    mesh: &HalfedgeMesh,
    b: usize,       // vertex B (placed at origin)
    h_out: usize,   // halfedge B → next (fan start)
    h_target: usize, // halfedge B → prev (fan end)
) -> Option<(Vec<usize>, Vec<[f64; 2]>)> {
    if mesh.origin(h_out) != b { return None; }
    if mesh.is_boundary_he(h_out) { return None; }

    let mut fan_he: Vec<usize> = vec![h_out];
    let mut verts_2d: Vec<[f64; 2]> = vec![[mesh.metric_edge_len(h_out), 0.0]];
    let mut cum = 0.0f64;
    let mut h = h_out;

    for _ in 0..512 {
        if mesh.is_boundary_he(h) { return None; }

        cum += mesh.metric_corner_angle(h);

        // CCW step around B: twin(prev(h))
        let ph = mesh.prev(h);
        if ph == INVALID { return None; }
        let hn = mesh.twin(ph);
        if hn == INVALID { return None; }

        let r = mesh.metric_edge_len(hn);
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
    b: usize,
    h_out: usize,
    h_target: usize,
) -> Option<(Vec<usize>, Vec<[f64; 2]>)> {
    if mesh.origin(h_out) != b { return None; }

    let mut fan_he: Vec<usize> = vec![h_out];
    let mut verts_2d: Vec<[f64; 2]> = vec![[mesh.metric_edge_len(h_out), 0.0]];
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

        cum -= mesh.metric_corner_angle(hn);

        let r = mesh.metric_edge_len(hn);
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
/// `a2` to `c2` with the **radial** fan edges.
///
/// Radial edges: from B=(0,0) to `verts_2d[k]` for k in `1..n-1`.
/// (k=0 and k=n-1 bound the fan — they are path edges, not interior.)
///
/// Each crossing at segment-parameter `t` along B→V_k becomes a
/// `PathPoint::Edge { v0: B, v1: V_k, t }`.  Multiple crossings are
/// returned sorted by `ray_s` (order A→C along the geodesic).
fn find_fan_crossings(
    mesh: &HalfedgeMesh,
    fan_he: &[usize],
    verts_2d: &[[f64; 2]],
    a2: [f64; 2],
    c2: [f64; 2],
) -> Vec<PathPoint> {
    let n = verts_2d.len();
    if n < 3 { return Vec::new(); }

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

/// For an Edge-point prev on edge (ev0, ev1), find the fan target halfedge
/// from vertex v. Returns the halfedge from v that should be the last in the
/// fan (the later one in CCW order within the triangle {v, ev0, ev1}).
fn fan_target_for_edge(mesh: &HalfedgeMesh, v: usize, ev0: usize, ev1: usize) -> Option<usize> {
    // Special case: edge includes the fan center vertex
    if ev0 == v {
        let h = find_halfedge(mesh, v, ev1);
        return if h != INVALID { Some(h) } else { None };
    }
    if ev1 == v {
        let h = find_halfedge(mesh, v, ev0);
        return if h != INVALID { Some(h) } else { None };
    }
    // Try v→ev0 in face containing ev1
    let h = find_halfedge(mesh, v, ev0);
    if h != INVALID && !mesh.is_boundary_he(h) && mesh.dest(mesh.next(h)) == ev1 {
        // Face is (v, ev0, ev1). CCW step from v→ev0 gives v→ev1.
        let hn = mesh.twin(mesh.prev(h));
        return if hn != INVALID { Some(hn) } else { None };
    }
    // Try v→ev1 in face containing ev0
    let h = find_halfedge(mesh, v, ev1);
    if h != INVALID && !mesh.is_boundary_he(h) && mesh.dest(mesh.next(h)) == ev0 {
        let hn = mesh.twin(mesh.prev(h));
        return if hn != INVALID { Some(hn) } else { None };
    }
    None
}

/// For an Edge-point next on edge (ev0, ev1), find the fan start halfedge
/// from vertex v. Returns the earlier halfedge from v in the triangle
/// {v, ev0, ev1} (the one directly in the face).
fn fan_start_for_edge(mesh: &HalfedgeMesh, v: usize, ev0: usize, ev1: usize) -> Option<usize> {
    // Special case: edge includes the fan center vertex
    if ev0 == v {
        let h = find_halfedge(mesh, v, ev1);
        return if h != INVALID { Some(h) } else { None };
    }
    if ev1 == v {
        let h = find_halfedge(mesh, v, ev0);
        return if h != INVALID { Some(h) } else { None };
    }
    let h = find_halfedge(mesh, v, ev0);
    if h != INVALID && !mesh.is_boundary_he(h) && mesh.dest(mesh.next(h)) == ev1 {
        return Some(h);
    }
    let h = find_halfedge(mesh, v, ev1);
    if h != INVALID && !mesh.is_boundary_he(h) && mesh.dest(mesh.next(h)) == ev0 {
        return Some(h);
    }
    None
}

/// Compute the 2-D position of an Edge point within the unfolded fan.
/// Finds ev0 and ev1 among the fan vertex destinations and interpolates.
/// `fan_center` is the vertex at the origin of the unfolded fan.
fn edge_point_2d_in_fan(
    fan_he: &[usize],
    verts_2d: &[[f64; 2]],
    mesh: &HalfedgeMesh,
    ev0: usize, ev1: usize, t: f64,
    fan_center: usize,
) -> Option<[f64; 2]> {
    // Special case: one edge vertex is the fan center (at origin)
    // Edge{v0, v1, t}: pos = pos(v0) + t*(pos(v1) - pos(v0))
    if ev0 == fan_center {
        // pos = origin + t*(pos(ev1) - origin) = t * pos(ev1)_2d
        let idx = fan_he.iter().position(|&h| mesh.dest(h) == ev1)?;
        let p = verts_2d[idx];
        return Some([t * p[0], t * p[1]]);
    }
    if ev1 == fan_center {
        // pos = pos(ev0) + t*(origin - pos(ev0)) = (1-t) * pos(ev0)_2d
        let idx = fan_he.iter().position(|&h| mesh.dest(h) == ev0)?;
        let p = verts_2d[idx];
        return Some([(1.0 - t) * p[0], (1.0 - t) * p[1]]);
    }
    // Normal case: both vertices in fan
    let mut idx0 = None;
    let mut idx1 = None;
    for (k, &h) in fan_he.iter().enumerate() {
        let d = mesh.dest(h);
        if d == ev0 { idx0 = Some(k); }
        if d == ev1 { idx1 = Some(k); }
    }
    let i0 = idx0?;
    let i1 = idx1?;
    let p0 = verts_2d[i0];
    let p1 = verts_2d[i1];
    Some([p0[0] + t * (p1[0] - p0[0]), p0[1] + t * (p1[1] - p0[1])])
}

/// Compute the 2-D position of a path endpoint in the unfolded fan.
fn endpoint_2d(
    pt: &PathPoint,
    fan_he: &[usize],
    verts_2d: &[[f64; 2]],
    mesh: &HalfedgeMesh,
    default_idx: usize,
    fan_center: usize,
) -> Option<[f64; 2]> {
    match pt {
        PathPoint::Vertex(_) => Some(verts_2d[default_idx]),
        PathPoint::Edge { v0, v1, t } => {
            edge_point_2d_in_fan(fan_he, verts_2d, mesh, *v0, *v1, *t, fan_center)
        }
    }
}

/// Try to straighten the path at Vertex `v` given its neighbours.
///
/// Returns the replacement PathPoints (0 = already straight, 1+ = edge crossings).
/// Handles both Vertex and Edge-point neighbours.
fn straighten_at_vertex(
    mesh: &HalfedgeMesh,
    prev_pt: &PathPoint,
    v: usize,
    next_pt: &PathPoint,
) -> Vec<PathPoint> {
    // Determine h_out (fan start) based on next_pt
    let h_out = match next_pt {
        PathPoint::Vertex(nv) => {
            if *nv == v { return Vec::new(); }
            let h = find_halfedge(mesh, v, *nv);
            if h == INVALID { return Vec::new(); }
            h
        }
        PathPoint::Edge { v0, v1, .. } => {
            match fan_start_for_edge(mesh, v, *v0, *v1) {
                Some(h) => h,
                None => return Vec::new(),
            }
        }
    };

    // Determine h_target (fan end) based on prev_pt
    let h_target = match prev_pt {
        PathPoint::Vertex(pv) => {
            if *pv == v { return Vec::new(); }
            let h_in = find_halfedge(mesh, *pv, v);
            if h_in == INVALID { return Vec::new(); }
            let ht = mesh.twin(h_in);
            if ht == INVALID { return Vec::new(); }
            ht
        }
        PathPoint::Edge { v0, v1, .. } => {
            match fan_target_for_edge(mesh, v, *v0, *v1) {
                Some(h) => h,
                None => return Vec::new(),
            }
        }
    };

    if h_out == h_target { return Vec::new(); }

    // Try CCW fan first (left turn at v)
    if let Some((fan_he, verts_2d)) = unfold_fan_ccw(mesh, v, h_out, h_target) {
        let n = verts_2d.len();
        if let (Some(a2), Some(c2)) = (
            endpoint_2d(prev_pt, &fan_he, &verts_2d, mesh, n - 1, v),
            endpoint_2d(next_pt, &fan_he, &verts_2d, mesh, 0, v),
        ) {
            let pts = find_fan_crossings(mesh, &fan_he, &verts_2d, a2, c2);
            if !pts.is_empty() { return pts; }
        }
    }

    // Try CW fan (right turn at v)
    if let Some((fan_he, verts_2d)) = unfold_fan_cw(mesh, v, h_out, h_target) {
        let n = verts_2d.len();
        if let (Some(a2), Some(c2)) = (
            endpoint_2d(prev_pt, &fan_he, &verts_2d, mesh, n - 1, v),
            endpoint_2d(next_pt, &fan_he, &verts_2d, mesh, 0, v),
        ) {
            let pts = find_fan_crossings(mesh, &fan_he, &verts_2d, a2, c2);
            if !pts.is_empty() { return pts; }
        }
    }

    Vec::new() // v is already locally straight (or a cone point)
}

// ────────────────────────────────────────────────────────────────
// Strip unfolding optimization
// ────────────────────────────────────────────────────────────────

/// Convert a vertex path into a sequence of crossed edges by collecting
/// the interior radial edges from each intermediate vertex's fan.
/// This allows strip unfolding to optimize paths that Dijkstra found
/// along mesh boundaries (where per-vertex shortening can't improve).
fn vertex_path_to_crossed_edges(
    mesh: &HalfedgeMesh,
    vpath: &[usize],
) -> Option<Vec<(usize, usize)>> {
    if vpath.len() < 3 { return None; }

    // For each intermediate vertex, compute BOTH fan directions and their
    // reversed interior edges. We'll pick the compatible direction later.
    let mut fan_options: Vec<(Vec<(usize, usize)>, Vec<(usize, usize)>)> = Vec::new();

    for i in 1..vpath.len() - 1 {
        let prev_v = vpath[i - 1];
        let v = vpath[i];
        let next_v = vpath[i + 1];
        if prev_v == v || next_v == v { return None; }

        let h_out = find_halfedge(mesh, v, next_v);
        if h_out == INVALID { return None; }
        let h_in = find_halfedge(mesh, prev_v, v);
        if h_in == INVALID { return None; }
        let h_target = mesh.twin(h_in);
        if h_target == INVALID { return None; }

        if h_out == h_target {
            fan_options.push((Vec::new(), Vec::new()));
            continue;
        }

        let extract = |fan_he: &[usize]| -> Vec<(usize, usize)> {
            (1..fan_he.len() - 1).rev()
                .map(|k| (mesh.origin(fan_he[k]), mesh.dest(fan_he[k])))
                .collect()
        };

        let ccw_edges = unfold_fan_ccw(mesh, v, h_out, h_target)
            .map(|(fan, _)| extract(&fan))
            .unwrap_or_default();
        let cw_edges = unfold_fan_cw(mesh, v, h_out, h_target)
            .map(|(fan, _)| extract(&fan))
            .unwrap_or_default();

        if ccw_edges.is_empty() && cw_edges.is_empty() {
            return None;
        }
        fan_options.push((ccw_edges, cw_edges));
    }

    // Build the crossed-edge sequence by choosing a compatible fan direction
    // at each step. Try both orderings for the first fan.
    let result = build_crossed_sequence(&fan_options, false)
        .or_else(|| build_crossed_sequence(&fan_options, true));

    result
}

/// Build a crossed-edge sequence from fan options, starting with CCW or CW.
fn build_crossed_sequence(
    fan_options: &[(Vec<(usize, usize)>, Vec<(usize, usize)>)],
    start_with_cw: bool,
) -> Option<Vec<(usize, usize)>> {
    let mut crossed = Vec::new();

    for (ccw, cw) in fan_options {
        if ccw.is_empty() && cw.is_empty() { continue; }

        if crossed.is_empty() {
            let (first, second) = if start_with_cw { (cw, ccw) } else { (ccw, cw) };
            if !first.is_empty() {
                crossed.extend_from_slice(first);
            } else if !second.is_empty() {
                crossed.extend_from_slice(second);
            }
        } else {
            let last = *crossed.last().unwrap();
            let cw_ok = !cw.is_empty() && shared_vertex(last, cw[0]).is_some();
            let ccw_ok = !ccw.is_empty() && shared_vertex(last, ccw[0]).is_some();
            if cw_ok {
                crossed.extend_from_slice(cw);
            } else if ccw_ok {
                crossed.extend_from_slice(ccw);
            } else {
                return None;
            }
        }
    }

    if crossed.is_empty() { return None; }

    // Verify consecutive edges share a vertex
    for w in crossed.windows(2) {
        if shared_vertex(w[0], w[1]).is_none() {
            return None;
        }
    }

    Some(crossed)
}

/// Find the vertex shared between two edges, if any.
fn shared_vertex(e1: (usize, usize), e2: (usize, usize)) -> Option<usize> {
    if e1.0 == e2.0 || e1.0 == e2.1 { Some(e1.0) }
    else if e1.1 == e2.0 || e1.1 == e2.1 { Some(e1.1) }
    else { None }
}

/// Place a new vertex across an edge in the 2D unfolded plane.
/// The new vertex is at distance `d_a` from `ea` and `d_b` from `eb`.
/// `above` selects which side of the edge to place it on.
fn place_vertex_2d(
    ea: [f64; 2], eb: [f64; 2],
    d_a: f64, d_b: f64,
    above: bool,
) -> [f64; 2] {
    let d_ab = ((eb[0] - ea[0]).powi(2) + (eb[1] - ea[1]).powi(2)).sqrt();
    if d_ab < 1e-15 { return ea; }
    let cos_a = (d_a * d_a + d_ab * d_ab - d_b * d_b) / (2.0 * d_a * d_ab);
    let a = cos_a.clamp(-1.0, 1.0).acos();
    let angle_ab = (eb[1] - ea[1]).atan2(eb[0] - ea[0]);
    let angle = if above { angle_ab + a } else { angle_ab - a };
    [ea[0] + d_a * angle.cos(), ea[1] + d_a * angle.sin()]
}

/// Optimize an open path segment by unfolding the entire triangle strip and
/// finding the straight-line geodesic. This gives the globally optimal edge
/// crossings for the given topological path (face sequence).
fn optimize_strip(mesh: &HalfedgeMesh, path: &GeodesicPath) -> GeodesicPath {
    if path.is_closed || path.points.len() < 3 {
        return path.clone();
    }

    let start_v = match path.points[0].as_vertex() {
        Some(v) => v,
        None => return path.clone(),
    };
    let end_v = match path.points.last().unwrap().as_vertex() {
        Some(v) => v,
        None => return path.clone(),
    };

    // Collect interior edge crossings
    let mut crossed: Vec<(usize, usize)> = Vec::new();
    for pt in &path.points[1..path.points.len() - 1] {
        match pt {
            PathPoint::Edge { v0, v1, .. } => crossed.push((*v0, *v1)),
            _ => return path.clone(),
        }
    }
    if crossed.is_empty() {
        return path.clone();
    }

    // Verify consecutive edges share a vertex (valid strip topology)
    for w in crossed.windows(2) {
        if shared_vertex(w[0], w[1]).is_none() {
            return path.clone();
        }
    }

    let nc = crossed.len();
    let mut pos_2d = std::collections::HashMap::<usize, [f64; 2]>::new();

    // Place first triangle: {start_v, crossed[0].0, crossed[0].1}
    let (e0a, e0b) = crossed[0];
    let d_s_a = metric_edge_len_between(mesh, start_v, e0a);
    let d_s_b = metric_edge_len_between(mesh, start_v, e0b);
    let d_a_b = metric_edge_len_between(mesh, e0a, e0b);
    if d_s_a.is_infinite() || d_s_b.is_infinite() || d_a_b.is_infinite() {
        return path.clone();
    }

    pos_2d.insert(start_v, [0.0, 0.0]);
    pos_2d.insert(e0a, [d_s_a, 0.0]);
    let cos_ang = (d_s_a * d_s_a + d_s_b * d_s_b - d_a_b * d_a_b) / (2.0 * d_s_a * d_s_b);
    let ang = cos_ang.clamp(-1.0, 1.0).acos();
    pos_2d.insert(e0b, [d_s_b * ang.cos(), d_s_b * ang.sin()]);

    let mut prev_v = start_v;

    for i in 0..nc {
        let (ev0, ev1) = crossed[i];

        let new_v = if i < nc - 1 {
            let (nv0, nv1) = crossed[i + 1];
            let shared = match shared_vertex((ev0, ev1), (nv0, nv1)) {
                Some(s) => s,
                None => return path.clone(),
            };
            if nv0 == shared { nv1 } else { nv0 }
        } else {
            end_v
        };

        if pos_2d.contains_key(&new_v) {
            if i < nc - 1 {
                let (nv0, nv1) = crossed[i + 1];
                let shared = shared_vertex((ev0, ev1), (nv0, nv1)).unwrap();
                prev_v = if ev0 == shared { ev1 } else { ev0 };
            }
            continue;
        }

        let ev0_2d = match pos_2d.get(&ev0) { Some(p) => *p, None => return path.clone() };
        let ev1_2d = match pos_2d.get(&ev1) { Some(p) => *p, None => return path.clone() };
        let prev_2d = match pos_2d.get(&prev_v) { Some(p) => *p, None => return path.clone() };

        let d0 = metric_edge_len_between(mesh, ev0, new_v);
        let d1 = metric_edge_len_between(mesh, ev1, new_v);
        if d0.is_infinite() || d1.is_infinite() { return path.clone(); }

        let edge_dir = sub2(ev1_2d, ev0_2d);
        let to_prev = sub2(prev_2d, ev0_2d);
        let prev_sign = cross2(edge_dir, to_prev);

        let try_above = prev_sign < 0.0;
        let placed = place_vertex_2d(ev0_2d, ev1_2d, d0, d1, try_above);
        let to_placed = sub2(placed, ev0_2d);
        let placed_sign = cross2(edge_dir, to_placed);

        let final_pos = if placed_sign * prev_sign > 0.0 {
            place_vertex_2d(ev0_2d, ev1_2d, d0, d1, !try_above)
        } else {
            placed
        };
        pos_2d.insert(new_v, final_pos);

        if i < nc - 1 {
            let (nv0, nv1) = crossed[i + 1];
            let shared = shared_vertex((ev0, ev1), (nv0, nv1)).unwrap();
            prev_v = if ev0 == shared { ev1 } else { ev0 };
        }
    }

    // Compute new crossings from the straight line
    let start_2d = pos_2d[&start_v];
    let end_2d = match pos_2d.get(&end_v) { Some(p) => *p, None => return path.clone() };

    let mut new_points = vec![PathPoint::Vertex(start_v)];
    for i in 0..nc {
        let (ev0, ev1) = crossed[i];
        let ev0_2d = pos_2d[&ev0];
        let ev1_2d = pos_2d[&ev1];

        if let Some((ray_s, edge_t)) = seg_seg_intersect(start_2d, end_2d, ev0_2d, ev1_2d) {
            let eps = 1e-12;
            if edge_t > eps && edge_t < 1.0 - eps && ray_s > eps && ray_s < 1.0 - eps {
                new_points.push(PathPoint::Edge { v0: ev0, v1: ev1, t: edge_t });
            }
        }
    }
    new_points.push(PathPoint::Vertex(end_v));

    let result = GeodesicPath { points: new_points, is_closed: false };

    // Only use the result if it's shorter
    let old_len = path.metric_length(mesh);
    let new_len = result.metric_length(mesh);
    if new_len < old_len - 1e-14 { result } else { path.clone() }
}

/// Apply strip unfolding in small overlapping windows along a vertex path.
/// Each window of W vertices is independently optimized. Small windows avoid
/// geometric distortion that occurs with long metric-weighted strips.
fn windowed_strip_optimize(mesh: &HalfedgeMesh, verts: &[usize]) -> Vec<PathPoint> {
    const W: usize = 8; // window size (vertices)
    let n = verts.len();
    if n < 3 {
        return verts.iter().map(|&v| PathPoint::Vertex(v)).collect();
    }

    // Try to optimize each window; collect (start_idx, optimized_segment) pairs
    let mut segments: Vec<(usize, Vec<PathPoint>)> = Vec::new();
    let mut i = 0;
    while i + 2 < n {
        let end = (i + W).min(n);
        let window = &verts[i..end];
        let mut improved = false;

        if window.len() >= 3 {
            if let Some(crossed) = vertex_path_to_crossed_edges(mesh, window) {
                if !crossed.is_empty() {
                    let ws = window[0];
                    let we = *window.last().unwrap();
                    let mut pp = vec![PathPoint::Vertex(ws)];
                    for &(v0, v1) in &crossed {
                        pp.push(PathPoint::Edge { v0, v1, t: 0.5 });
                    }
                    pp.push(PathPoint::Vertex(we));
                    let proxy = GeodesicPath { points: pp, is_closed: false };
                    let opt = optimize_strip(mesh, &proxy);
                    let orig = GeodesicPath {
                        points: window.iter().map(|&v| PathPoint::Vertex(v)).collect(),
                        is_closed: false,
                    };
                    if opt.metric_length(mesh) < orig.metric_length(mesh) - 1e-14 {
                        segments.push((i, opt.points));
                        i = end - 1; // next window starts at the last vertex of this one
                        improved = true;
                    }
                }
            }
        }
        if !improved {
            i += 1;
        }
    }

    // Merge: build the output path from original vertices + optimized segments
    if segments.is_empty() {
        return verts.iter().map(|&v| PathPoint::Vertex(v)).collect();
    }

    let mut result: Vec<PathPoint> = Vec::new();
    let mut cursor = 0usize;
    for (seg_start, seg_pts) in &segments {
        // Add original vertices from cursor to seg_start
        for j in cursor..*seg_start {
            result.push(PathPoint::Vertex(verts[j]));
        }
        // Add the optimized segment (including its start vertex)
        result.extend_from_slice(seg_pts);
        // The segment ends at verts[seg_start + W - 1] or similar;
        // find the last vertex of the segment to set cursor
        if let Some(last_v) = seg_pts.last().and_then(|p| p.as_vertex()) {
            cursor = verts.iter().position(|&v| v == last_v).unwrap_or(n);
        } else {
            cursor = n; // shouldn't happen
        }
    }
    // Add remaining vertices after the last segment
    for j in cursor..n {
        // Avoid duplicating the last vertex of the last segment
        if !result.is_empty() {
            if let Some(PathPoint::Vertex(last_v)) = result.last() {
                if *last_v == verts[j] { continue; }
            }
        }
        result.push(PathPoint::Vertex(verts[j]));
    }

    if result.len() >= 2 { result } else { verts.iter().map(|&v| PathPoint::Vertex(v)).collect() }
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
///
/// When the mesh is not intrinsically Delaunay, non-Delaunay edges can trap
/// the shortening in local optima. We first flip edges on a **copy** of the
/// mesh to make it Delaunay, run the shortening on the Delaunay copy, then
/// translate the result back to the original mesh's edge parameterization.
pub fn shorten_path(
    mesh: &HalfedgeMesh,
    path: &GeodesicPath,
    config: &FlipConfig,
) -> GeodesicPath {
    // For Euclidean metric: flip to intrinsic Delaunay and shorten on the
    // flipped copy. Intrinsic flipping is exact for Euclidean geometry.
    // For non-Euclidean metrics: skip flipping (the metric-weighted strip
    // unfolding and metric_length approximation cause accuracy issues).
    if mesh.is_euclidean() {
        let mut dm = mesh.clone();
        let n_flips = dm.flip_to_delaunay(dm.n_edges() * 4);
        if n_flips > 0 {
            let result_on_dm = shorten_path_inner(&dm, path, config);
            let translated = translate_path_to_original(mesh, &dm, &result_on_dm);
            return translated;
        }
    }
    shorten_path_inner(mesh, path, config)
}

/// Translate a path from a Delaunay-flipped mesh back to the original mesh.
/// Edge crossings on flipped edges are converted to crossings on the
/// original mesh edges by finding the 3-D position and locating the
/// enclosing original triangle.
fn translate_path_to_original(
    orig: &HalfedgeMesh,
    flipped: &HalfedgeMesh,
    path: &GeodesicPath,
) -> GeodesicPath {
    let mut new_points = Vec::with_capacity(path.points.len());
    for pt in &path.points {
        match pt {
            PathPoint::Vertex(v) => {
                new_points.push(PathPoint::Vertex(*v));
            }
            PathPoint::Edge { v0, v1, t } => {
                // Check if this edge exists in the original mesh
                let h = find_halfedge(orig, *v0, *v1);
                if h != INVALID {
                    // Edge exists in original — keep as-is
                    new_points.push(pt.clone());
                } else {
                    // Edge was created by flipping — find 3-D position and
                    // locate the enclosing edge in the original mesh.
                    let pos = pt.position(flipped);
                    if let Some(new_pt) = locate_point_on_original(orig, &pos, *v0, *v1) {
                        new_points.push(new_pt);
                    } else {
                        // Fallback: use nearest vertex
                        let dv0 = crate::mesh::dist3(&pos, &orig.position(*v0));
                        let dv1 = crate::mesh::dist3(&pos, &orig.position(*v1));
                        new_points.push(PathPoint::Vertex(if dv0 <= dv1 { *v0 } else { *v1 }));
                    }
                }
            }
        }
    }
    GeodesicPath { points: new_points, is_closed: path.is_closed }
}

/// Find which edge/face of the original mesh contains a 3-D point that lies
/// on a flipped edge (v0, v1). The flipped edge (v0,v1) replaced the
/// original edge between the two triangles sharing v0 and v1's quad.
fn locate_point_on_original(
    orig: &HalfedgeMesh,
    pos: &[f64; 3],
    v0: usize,
    v1: usize,
) -> Option<PathPoint> {
    // The flipped edge (v0, v1) was the diagonal of a quad. In the original
    // mesh, this quad had the OTHER diagonal. Find it: look for a vertex
    // adjacent to both v0 and v1 in the original mesh.
    // These shared neighbors are the other two vertices of the original quad.
    let mut shared = Vec::new();
    for h in orig.outgoing_halfedges(v0) {
        let w = orig.dest(h);
        if w == v1 { continue; }
        for h2 in orig.outgoing_halfedges(v1) {
            if orig.dest(h2) == w {
                shared.push(w);
                break;
            }
        }
    }

    // The original edge was between two shared neighbors.
    // The flipped edge (v0,v1) crosses the original edge (shared[0], shared[1]).
    if shared.len() >= 2 {
        let a = shared[0];
        let b = shared[1];
        let h = find_halfedge(orig, a, b);
        if h != INVALID {
            // Project pos onto edge a→b
            let pa = orig.position(a);
            let pb = orig.position(b);
            let ab = crate::mesh::sub3(&pb, &pa);
            let ap = crate::mesh::sub3(pos, &pa);
            let len_sq = crate::mesh::dot3(&ab, &ab);
            if len_sq > 1e-30 {
                let t = (crate::mesh::dot3(&ap, &ab) / len_sq).clamp(0.0, 1.0);
                return Some(PathPoint::Edge { v0: a, v1: b, t });
            }
        }
    }
    None
}

/// Inner shortening logic (operates on a single mesh without flipping).
fn shorten_path_inner(
    mesh: &HalfedgeMesh,
    path: &GeodesicPath,
    config: &FlipConfig,
) -> GeodesicPath {
    let is_closed = path.is_closed;

    if path.points.len() < 3 {
        return path.clone();
    }

    // Pre-optimization: for vertex paths, apply strip unfolding to convert
    // boundary vertices to edge crossings.
    // - Euclidean: full-path strip (exact, handles all boundary vertices)
    // - Non-Euclidean: windowed strip (limits metric distortion)
    let mut pts = path.points.clone();
    if !is_closed {
        let verts: Vec<usize> = pts.iter().filter_map(|p| p.as_vertex()).collect();
        if verts.len() == pts.len() && verts.len() >= 3 {
            if mesh.is_euclidean() {
                // Full-path strip unfolding (exact for Euclidean)
                if let Some(crossed) = vertex_path_to_crossed_edges(mesh, &verts) {
                    if !crossed.is_empty() {
                        let start_v = verts[0];
                        let end_v = *verts.last().unwrap();
                        let mut pp = vec![PathPoint::Vertex(start_v)];
                        for &(v0, v1) in &crossed {
                            pp.push(PathPoint::Edge { v0, v1, t: 0.5 });
                        }
                        pp.push(PathPoint::Vertex(end_v));
                        let proxy = GeodesicPath { points: pp, is_closed: false };
                        let opt = optimize_strip(mesh, &proxy);
                        if opt.metric_length(mesh) < path.metric_length(mesh) - 1e-14 {
                            pts = opt.points;
                        }
                    }
                }
            }
        }
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
            .metric_length(mesh);

        if cur_len < 1e-15 { break; }
        let rel_dec = (prev_len - cur_len) / prev_len.max(cur_len);
        if !changed || rel_dec.abs() < config.rel_tol { break; }
        prev_len = cur_len;
    }

    let mut result = GeodesicPath { points: pts, is_closed };

    // Post-process: strip unfolding for all-edge-crossing paths
    if !is_closed {
        result = optimize_strip(mesh, &result);
    }

    // Post-process: windowed strip unfolding for remaining vertex runs.
    // After per-vertex shortening, some vertices (especially on mesh boundaries)
    // remain unshortened. Find contiguous runs of Vertex points and optimize them.
    if !is_closed {
        let pts = &result.points;
        // Find runs of consecutive Vertex points of length >= 3
        let mut new_pts: Vec<PathPoint> = Vec::new();
        let mut run_start = None;
        for i in 0..=pts.len() {
            let is_vertex = i < pts.len() && pts[i].as_vertex().is_some();
            if is_vertex {
                if run_start.is_none() { run_start = Some(i); }
            } else {
                if let Some(start) = run_start {
                    let run = &pts[start..i];
                    let run_verts: Vec<usize> = run.iter()
                        .filter_map(|p| p.as_vertex()).collect();
                    if run_verts.len() >= 3 {
                        let optimized = windowed_strip_optimize(mesh, &run_verts);
                        // Skip first if it duplicates the last non-vertex point
                        if start > 0 && !new_pts.is_empty() {
                            new_pts.extend_from_slice(&optimized);
                        } else {
                            new_pts.extend(optimized);
                        }
                    } else {
                        new_pts.extend_from_slice(run);
                    }
                    run_start = None;
                }
                if i < pts.len() {
                    new_pts.push(pts[i].clone());
                }
            }
        }
        if new_pts.len() >= 2 {
            result = GeodesicPath { points: new_pts, is_closed: false };
        }
    }

    result
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
    // Angles at the OPPOSITE vertices (not the edge endpoints)
    let alpha = mesh.corner_angle(mesh.next(mesh.next(h)));
    let beta  = mesh.corner_angle(mesh.next(mesh.next(t)));
    alpha + beta <= std::f64::consts::PI + 1e-10
}
