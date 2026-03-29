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
    /// For Euclidean metric, uses exact edge/segment lengths.
    /// For non-Euclidean metrics, integrates the speed field along each segment
    /// via barycentric interpolation. This gives consistent results regardless
    /// of whether the path has edge crossings or only vertex points.
    pub fn metric_length(&self, mesh: &HalfedgeMesh) -> f64 {
        if mesh.is_euclidean() {
            return self.euclidean_length(mesh);
        }
        let positions: Vec<[f64; 3]> = self.points.iter()
            .map(|p| p.position(mesh)).collect();
        let n = positions.len();
        if n < 2 { return 0.0; }
        let seg_count = if self.is_closed { n } else { n - 1 };
        let mut total = 0.0;
        for s in 0..seg_count {
            let j = (s + 1) % n;
            total += integrate_metric_segment(mesh, &positions[s], &positions[j]);
        }
        total
    }
}

// ────────────────────────────────────────────────────────────────
// Internal mesh helpers
// ────────────────────────────────────────────────────────────────

/// Integrate the metric length along a straight-line 3-D segment from `p` to `q`.
/// Subdivides the segment, locates each sample in a mesh triangle, and
/// interpolates the per-vertex metric ratio barycentrically.
fn integrate_metric_segment(mesh: &HalfedgeMesh, p: &[f64; 3], q: &[f64; 3]) -> f64 {
    let euc = dist3(p, q);
    if euc < 1e-30 { return 0.0; }

    let avg_edge = if mesh.n_edges() > 0 {
        mesh.edge_lengths.iter().sum::<f64>() / mesh.n_edges() as f64
    } else {
        euc
    };
    let n_steps = ((euc / avg_edge).ceil() as usize).max(1).min(500);
    let dt = 1.0 / n_steps as f64;
    let step_len = euc * dt;

    // Precompute per-vertex slowness (1/speed) for barycentric interpolation
    let slowness = vertex_slowness_array(mesh);

    // Seed face search from a vertex near the start
    let start_v = nearest_vertex_brute(mesh, p);
    let mut hint_face = mesh.vert_he[start_v];
    if hint_face != INVALID { hint_face = mesh.he_face[hint_face]; }

    let mut total = 0.0;
    for i in 0..n_steps {
        let t = (i as f64 + 0.5) * dt;
        let pt = [
            p[0] + t * (q[0] - p[0]),
            p[1] + t * (q[1] - p[1]),
            p[2] + t * (q[2] - p[2]),
        ];
        let (ratio, face) = barycentric_metric_ratio(mesh, &pt, &slowness, hint_face);
        if face != INVALID { hint_face = face; }
        total += step_len * ratio;
    }
    total
}

/// Build per-vertex slowness array (metric/euclidean ratio ≈ 1/speed).
fn vertex_slowness_array(mesh: &HalfedgeMesh) -> Vec<f64> {
    match &mesh.metric {
        crate::mesh::Metric::Euclidean => vec![1.0; mesh.n_verts()],
        crate::mesh::Metric::IsotropicSpeed(speeds) => {
            speeds.iter().map(|&s| if s > 1e-30 { 1.0 / s } else { 1e30 }).collect()
        }
        crate::mesh::Metric::FiberTensor { speeds_transverse, .. } => {
            // For scalar ratio approximation, use transverse speed
            speeds_transverse.iter().map(|&s| if s > 1e-30 { 1.0 / s } else { 1e30 }).collect()
        }
    }
}

/// Locate the face containing `pt` (walking from `hint`), compute the
/// barycentric interpolation of the slowness field, and return (ratio, face).
fn barycentric_metric_ratio(
    mesh: &HalfedgeMesh,
    pt: &[f64; 3],
    slowness: &[f64],
    hint: usize,
) -> (f64, usize) {
    let face = locate_face(mesh, pt, hint);
    if face == INVALID {
        // Fallback: nearest vertex
        let v = nearest_vertex_brute(mesh, pt);
        return (slowness[v], INVALID);
    }
    let h0 = mesh.face_he[face];
    let h1 = mesh.he_next[h0];
    let h2 = mesh.he_next[h1];
    let va = mesh.he_origin[h0];
    let vb = mesh.he_origin[h1];
    let vc = mesh.he_origin[h2];
    let pa = mesh.position(va);
    let pb = mesh.position(vb);
    let pc = mesh.position(vc);

    let (u, v, w) = barycentric_coords(pt, &pa, &pb, &pc);
    let ratio = u * slowness[va] + v * slowness[vb] + w * slowness[vc];
    (ratio.max(0.0), face)
}

/// Barycentric coordinates of point `p` in triangle (a, b, c).
/// Returns (u, v, w) with u+v+w ≈ 1, clamped to [0,1].
fn barycentric_coords(
    p: &[f64; 3], a: &[f64; 3], b: &[f64; 3], c: &[f64; 3],
) -> (f64, f64, f64) {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let ap = sub3(p, a);
    let d00 = dot3(&ab, &ab);
    let d01 = dot3(&ab, &ac);
    let d11 = dot3(&ac, &ac);
    let d20 = dot3(&ap, &ab);
    let d21 = dot3(&ap, &ac);
    let denom = d00 * d11 - d01 * d01;
    if denom.abs() < 1e-30 {
        return (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    }
    let v = (d11 * d20 - d01 * d21) / denom;
    let w = (d00 * d21 - d01 * d20) / denom;
    let u = 1.0 - v - w;
    // Clamp negatives and renormalize so u+v+w = 1
    let uc = u.max(0.0);
    let vc = v.max(0.0);
    let wc = w.max(0.0);
    let s = uc + vc + wc;
    if s > 1e-30 { (uc / s, vc / s, wc / s) } else { (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0) }
}

/// Walk across faces to locate the face containing `pt`, starting from `hint`.
fn locate_face(mesh: &HalfedgeMesh, pt: &[f64; 3], hint: usize) -> usize {
    if hint >= mesh.n_faces() {
        // No valid hint — find any face from nearest vertex
        let v = nearest_vertex_brute(mesh, pt);
        for h in mesh.outgoing_halfedges(v) {
            let f = mesh.he_face[h];
            if f != INVALID { return locate_face(mesh, pt, f); }
        }
        return INVALID;
    }

    let mut face = hint;
    for _ in 0..mesh.n_faces().min(200) {
        let h0 = mesh.face_he[face];
        let h1 = mesh.he_next[h0];
        let h2 = mesh.he_next[h1];

        let pa = mesh.position(mesh.he_origin[h0]);
        let pb = mesh.position(mesh.he_origin[h1]);
        let pc = mesh.position(mesh.he_origin[h2]);

        let (u, v, w) = barycentric_coords(pt, &pa, &pb, &pc);

        // If inside (or close enough), return this face
        let eps = -1e-6;
        if u >= eps && v >= eps && w >= eps {
            return face;
        }

        // Walk toward the point: cross the edge opposite the most negative coord
        let cross_h = if u < v && u < w { h1 }      // opposite to a → cross b-c
                      else if v < w     { h2 }      // opposite to b → cross c-a
                      else              { h0 };      // opposite to c → cross a-b
        let twin = mesh.he_twin[cross_h];
        let next_face = mesh.he_face[twin];
        if next_face == INVALID { return face; } // hit boundary, use current
        face = next_face;
    }
    face
}

/// Brute-force nearest vertex (used to seed face walks).
fn nearest_vertex_brute(mesh: &HalfedgeMesh, pos: &[f64; 3]) -> usize {
    let mut best_d2 = f64::INFINITY;
    let mut best_v = 0usize;
    for v in 0..mesh.n_verts() {
        let vp = mesh.position(v);
        let d2 = (pos[0]-vp[0]).powi(2) + (pos[1]-vp[1]).powi(2) + (pos[2]-vp[2]).powi(2);
        if d2 < best_d2 { best_d2 = d2; best_v = v; }
    }
    best_v
}

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

/// Shorten a geodesic path using FlipOut (Sharp & Crane 2020) on an
/// intrinsic triangulation, then extract edge crossings on the original mesh.
/// Falls back to Polthier-Schmies if FlipOut can't improve the path.
pub fn shorten_path(
    mesh: &HalfedgeMesh,
    path: &GeodesicPath,
    config: &FlipConfig,
) -> GeodesicPath {
    // Step 1: Polthier-Schmies shortening (fast, handles edge crossings)
    let ps_result = shorten_path_inner(mesh, path, config);

    // Step 2: FlipOut on the Dijkstra vertex path (handles non-straight
    // vertices that Polthier-Schmies can't improve, e.g. path turns).
    // Skip for Euclidean: PS + strip unfolding is already exact.
    if mesh.is_euclidean() { return ps_result; }
    if path.points.len() < 3 { return ps_result; }
    let verts: Vec<usize> = path.points.iter().filter_map(|p| p.as_vertex()).collect();
    if verts.len() != path.points.len() { return ps_result; }

    let mut imesh = mesh.clone();
    if !mesh.is_euclidean() {
        imesh.init_metric_edge_lengths();
    }

    let mut edge_path: Vec<usize> = Vec::new();
    for w in verts.windows(2) {
        let h = find_halfedge(&imesh, w[0], w[1]);
        if h == INVALID { return ps_result; }
        edge_path.push(h);
    }

    let initial_len: f64 = edge_path.iter().map(|&h| imesh.edge_len(h)).sum();
    let mut failed_vertices = std::collections::HashSet::new();
    let max_iters = (edge_path.len() * 10).min(config.max_iterations);
    for _iter in 0..max_iters {
        let n_edges = edge_path.len();
        if n_edges < 2 { break; }

        let mut worst_idx = None;
        let mut worst_angle = std::f64::consts::PI - 1e-10;
        for i in 0..n_edges - 1 {
            let b = imesh.dest(edge_path[i]);
            if failed_vertices.contains(&b) { continue; }
            let angle = smaller_wedge_angle(&imesh, b, edge_path[i], edge_path[i + 1]);
            if angle < worst_angle {
                worst_angle = angle;
                worst_idx = Some(i);
            }
        }
        let idx = match worst_idx { Some(i) => i, None => break };

        let h_in = edge_path[idx];
        let h_out = if idx + 1 < n_edges { edge_path[idx + 1] } else { edge_path[0] };
        let b = imesh.dest(h_in);

        if let Some(outer_arc) = flip_out_at_vertex(&mut imesh, b, h_in, h_out) {
            let arc_len: f64 = outer_arc.iter().map(|&h| imesh.edge_len(h)).sum();
            let replaced_len = imesh.edge_len(h_in) + imesh.edge_len(h_out);
            if arc_len >= replaced_len {
                // Outer arc is NOT shorter — skip this vertex
                failed_vertices.insert(b);
                continue;
            }
            if idx + 1 < n_edges {
                edge_path.splice(idx..=idx + 1, outer_arc);
            } else {
                edge_path.pop();
                edge_path.splice(0..1, outer_arc);
            }
            failed_vertices.clear();
        } else {
            failed_vertices.insert(b);
        }
    }

    let final_len: f64 = edge_path.iter().map(|&h| imesh.edge_len(h)).sum();
    if final_len < initial_len - 1e-14 {
        // FlipOut improved the path: extract and compare with PS result
        let flipout_result = extract_edge_path_on_original(mesh, &imesh, &edge_path, path.is_closed);
        let fo_len = flipout_result.metric_length(mesh);
        let ps_len = ps_result.metric_length(mesh);
        if fo_len < ps_len { flipout_result } else { ps_result }
    } else {
        ps_result
    }
}


fn smaller_wedge_angle(mesh: &HalfedgeMesh, b: usize, h_in: usize, h_out: usize) -> f64 {
    // Walk CCW from h_out to twin(h_in) = b→origin(h_in)
    let h_target = mesh.twin(h_in); // halfedge from b toward prev vertex
    let mut angle = 0.0f64;
    let mut h = h_out;
    for _ in 0..512 {
        if mesh.is_boundary_he(h) { return f64::INFINITY; }
        angle += mesh.corner_angle(h);
        let ph = mesh.prev(h);
        if ph == INVALID { return f64::INFINITY; }
        let hn = mesh.twin(ph);
        if hn == INVALID { return f64::INFINITY; }
        if hn == h_target { break; }
        h = hn;
        if angle > 2.0 * std::f64::consts::PI { return f64::INFINITY; }
    }
    // Also compute CW angle
    let total = mesh.vertex_angle_sum(b);
    let other = total - angle;
    angle.min(other)
}

/// FlipOut at vertex b (Algorithm 1 from Sharp & Crane 2020).
///
/// Repeatedly flips edges from b into the smaller wedge until all
/// outer arc angles β_i ≥ π, then returns the outer arc as halfedges.
fn flip_out_at_vertex(
    mesh: &mut HalfedgeMesh,
    b: usize,
    h_in: usize,   // halfedge prev → b
    h_out: usize,   // halfedge b → next
) -> Option<Vec<usize>> {
    let h_target = mesh.twin(h_in); // b → prev

    // Iteratively flip edges out of the wedge
    for _round in 0..500 {
        let wedge = collect_wedge_ccw(mesh, b, h_out, h_target)?;
        let k = wedge.len() - 1; // degree = number of triangles
        if k <= 1 { break; }     // single triangle or empty → done

        // Check each interior vertex n_j (wedge[1]..wedge[k-1]) for β < π.
        // β < π in the unfolded plane ⟺ sum of angles at n_j in the two
        // adjacent wedge faces > π (the edge b→n_j is "non-Delaunay" in the wedge).
        let mut did_flip = false;
        for j in 1..k {
            let alpha_left = mesh.corner_angle(mesh.prev(wedge[j - 1]));
            let alpha_right = mesh.corner_angle(mesh.next(wedge[j]));
            if alpha_left + alpha_right > std::f64::consts::PI + 1e-10 {
                // β_j < π: flip edge b→n_j out of the wedge
                // After flip, b→n_j becomes n_{j-1}→n_{j+1}. This removes
                // n_j from the wedge and decreases its degree by 1.
                // h_out and h_target must remain valid after the flip.
                let he = wedge[j];
                // Check that flipping won't invalidate h_out or h_target
                // (they could be in adjacent faces that get modified)
                if mesh.flip_edge(he) {
                    // Verify h_out and h_target are still valid
                    if mesh.origin(h_out) != b || mesh.origin(h_target) != b {
                        // Flip broke our references — abort
                        return None;
                    }
                    did_flip = true;
                    break;
                }
            }
        }
        if !did_flip { break; } // all β ≥ π → wedge is convex
    }

    // Collect the outer arc: for each face in the final wedge,
    // the outer edge is next(wedge[j]) = edge from n_j to n_{j+1}.
    let wedge = collect_wedge_ccw(mesh, b, h_out, h_target)?;
    let k = wedge.len() - 1;
    let mut outer_arc = Vec::new();
    for j in 0..k {
        outer_arc.push(mesh.next(wedge[j]));
    }
    // outer_arc goes: c→n1, n1→n2, ..., n_{k-1}→a
    // But we need a→...→c. Reverse and use twins.
    outer_arc.reverse();
    outer_arc = outer_arc.iter().map(|&h| mesh.twin(h)).collect();

    // Verify connectivity: each edge must connect to the next
    let a = mesh.origin(h_in);
    let c = mesh.dest(h_out);
    if outer_arc.is_empty() { return None; }
    if mesh.origin(outer_arc[0]) != a { return None; }
    if mesh.dest(*outer_arc.last().unwrap()) != c { return None; }
    for w in outer_arc.windows(2) {
        if mesh.dest(w[0]) != mesh.origin(w[1]) { return None; }
    }

    Some(outer_arc)
}

/// Collect CCW wedge halfedges from b (halfedges from b to outer vertices).
fn collect_wedge_ccw(
    mesh: &HalfedgeMesh, b: usize, h_out: usize, h_target: usize,
) -> Option<Vec<usize>> {
    let mut result = vec![h_out];
    let mut h = h_out;
    for _ in 0..512 {
        if mesh.is_boundary_he(h) { return None; }
        let ph = mesh.prev(h);
        if ph == INVALID { return None; }
        let hn = mesh.twin(ph);
        if hn == INVALID { return None; }
        result.push(hn);
        if hn == h_target { return Some(result); }
        h = hn;
    }
    None
}

/// Collect CW wedge halfedges from b.
fn collect_wedge_cw(
    mesh: &HalfedgeMesh, b: usize, h_out: usize, h_target: usize,
) -> Option<Vec<usize>> {
    let mut result = vec![h_out];
    let mut h = h_out;
    for _ in 0..512 {
        let th = mesh.twin(h);
        if th == INVALID { return None; }
        if mesh.is_boundary_he(th) { return None; }
        let hn = mesh.next(th);
        if hn == INVALID { return None; }
        if mesh.origin(hn) != b { return None; }
        result.push(hn);
        if hn == h_target { return Some(result); }
        h = hn;
    }
    None
}

/// Extract an intrinsic edge path as PathPoints on the original mesh.
/// Each edge in the intrinsic mesh either exists in the original (keep as Vertex
/// endpoints) or was created by flipping (convert to Edge crossing).
fn extract_edge_path_on_original(
    orig: &HalfedgeMesh,
    imesh: &HalfedgeMesh,
    edge_path: &[usize],
    is_closed: bool,
) -> GeodesicPath {
    if edge_path.is_empty() {
        return GeodesicPath { points: Vec::new(), is_closed };
    }

    let mut points = Vec::new();
    // Add the first vertex
    points.push(PathPoint::Vertex(imesh.origin(edge_path[0])));

    for &he in edge_path {
        let v0 = imesh.origin(he);
        let v1 = imesh.dest(he);

        // Check if this edge exists in the original mesh
        let h_orig = find_halfedge(orig, v0, v1);
        if h_orig != INVALID {
            // Edge exists — just add the destination vertex
            points.push(PathPoint::Vertex(v1));
        } else {
            // Edge was created by flipping — find the 3D position of the midpoint
            // and project onto the original mesh
            let pos = imesh.position(v0); // approximate: use v0 position
            let pos2 = imesh.position(v1);
            let mid = [
                (pos[0] + pos2[0]) * 0.5,
                (pos[1] + pos2[1]) * 0.5,
                (pos[2] + pos2[2]) * 0.5,
            ];
            if let Some(pt) = locate_point_on_original(orig, &mid, v0, v1) {
                points.push(pt);
            }
            points.push(PathPoint::Vertex(v1));
        }
    }

    // For closed paths, remove the duplicate last vertex
    if is_closed && points.len() >= 2 {
        if let (Some(PathPoint::Vertex(f)), Some(PathPoint::Vertex(l))) =
            (points.first(), points.last()) {
            if f == l { points.pop(); }
        }
    }

    GeodesicPath { points, is_closed }
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

    // Post-process: for each contiguous run of ≥3 Vertex points, apply
    // strip unfolding. Include the adjacent Edge-crossing points so the
    // optimizer sees where the path comes from / goes to (otherwise collinear
    // boundary runs can't be straightened).
    if !is_closed && result.points.len() >= 3 {
        let pts = result.points.clone();
        let n = pts.len();

        // Identify vertex runs (start_index, length)
        let mut runs: Vec<(usize, usize)> = Vec::new();
        let mut run_start: Option<usize> = None;
        for i in 0..=n {
            let is_v = i < n && pts[i].as_vertex().is_some();
            if is_v {
                if run_start.is_none() { run_start = Some(i); }
            } else if let Some(s) = run_start {
                if i - s >= 3 { runs.push((s, i - s)); }
                run_start = None;
            }
        }

        // Process runs in reverse order (so splice indices stay valid)
        for &(run_start_idx, run_len) in runs.iter().rev() {
            let run_end_idx = run_start_idx + run_len; // exclusive
            let run_verts: Vec<usize> = pts[run_start_idx..run_end_idx]
                .iter().filter_map(|p| p.as_vertex()).collect();
            if run_verts.len() < 3 { continue; }

            // Build extended vertex list: the run's vertices plus one
            // adjacent vertex on each side extracted from neighboring
            // edge crossings (their v0 gives the hub vertex).
            let mut ext_verts = run_verts.clone();
            // Prepend: scan backward to find a non-collinear hub vertex
            {
                let first_v = ext_verts[0];
                let mut found = false;
                for k in (run_start_idx.saturating_sub(6)..run_start_idx).rev() {
                    match &pts[k] {
                        PathPoint::Edge { v0, v1, .. } => {
                            for &candidate in &[*v1, *v0] {
                                if candidate != first_v {
                                    let h = find_halfedge(mesh, candidate, first_v);
                                    if h != INVALID {
                                        ext_verts.insert(0, candidate);
                                        found = true;
                                        break;
                                    }
                                }
                            }
                            if found { break; }
                        }
                        PathPoint::Vertex(v) => {
                            if *v != first_v {
                                ext_verts.insert(0, *v);
                                found = true;
                            }
                            break;
                        }
                    }
                }
            }
            // Append: scan forward past any edge crossings to find a
            // vertex that gives directional context (preferably non-collinear
            // with the boundary run). Try v0, v1 of each edge crossing,
            // then the next Vertex point.
            {
                let last_v = *ext_verts.last().unwrap();
                let mut found = false;
                for k in run_end_idx..n.min(run_end_idx + 6) {
                    match &pts[k] {
                        PathPoint::Edge { v0, v1, .. } => {
                            // Prefer the vertex that's NOT on the boundary
                            // (gives diagonal direction info)
                            for &candidate in &[*v1, *v0] {
                                if candidate != last_v {
                                    let h = find_halfedge(mesh, last_v, candidate);
                                    if h != INVALID {
                                        ext_verts.push(candidate);
                                        found = true;
                                        break;
                                    }
                                }
                            }
                            if found { break; }
                        }
                        PathPoint::Vertex(v) => {
                            if *v != last_v {
                                ext_verts.push(*v);
                                found = true;
                            }
                            break;
                        }
                    }
                }
            }
            if ext_verts.len() < 3 { continue; }
            // De-dup consecutive identical vertices
            ext_verts.dedup();
            if ext_verts.len() < 3 { continue; }
            let ext_start = if run_start_idx > 0 { run_start_idx - 1 } else { run_start_idx };
            let ext_end = if run_end_idx < n { run_end_idx + 1 } else { run_end_idx };

            // Build crossed edges and optimize the vertex sub-path
            if let Some(crossed) = vertex_path_to_crossed_edges(mesh, &ext_verts) {
                if !crossed.is_empty() {
                    let ws = ext_verts[0];
                    let we = *ext_verts.last().unwrap();
                    let mut pp = vec![PathPoint::Vertex(ws)];
                    for &(v0, v1) in &crossed {
                        pp.push(PathPoint::Edge { v0, v1, t: 0.5 });
                    }
                    pp.push(PathPoint::Vertex(we));
                    let proxy = GeodesicPath { points: pp, is_closed: false };
                    // Use Euclidean geometry for the strip optimization:
                    // boundary vertex runs are in fast zones where metric ≈ Euclidean.
                    // Metric strip distortion would reject valid improvements.
                    let mut euc_mesh = mesh.clone();
                    euc_mesh.clear_metric();
                    let opt = optimize_strip(&euc_mesh, &proxy);
                    let orig_seg = GeodesicPath {
                        points: pts[ext_start..ext_end].to_vec(),
                        is_closed: false,
                    };
                    // Compare using metric_length on the ORIGINAL mesh
                    // (the strip is Euclidean-correct, but the acceptance
                    // criterion must respect the speed field)
                    if opt.metric_length(mesh) < orig_seg.metric_length(mesh) - 1e-14 {
                        let mut new_pts = result.points[..ext_start].to_vec();
                        new_pts.extend(opt.points);
                        new_pts.extend_from_slice(&result.points[ext_end..]);
                        result = GeodesicPath { points: new_pts, is_closed: false };
                    }
                }
            }
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
