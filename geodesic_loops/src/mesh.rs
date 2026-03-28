/// Half-edge based triangulated surface mesh.
///
/// Halfedges are paired: halfedge `h` and `he_twin[h]` share an edge.
/// The face to the left of `h` is `he_face[h]` (INVALID for boundary halfedges).
/// `he_next[h]` is the next halfedge around the same face (CCW order).
/// `he_origin[h]` is the source vertex of `h`.

pub const INVALID: usize = usize::MAX;

/// Metric mode for geodesic computation.
#[derive(Debug, Clone)]
pub enum Metric {
    /// Standard Euclidean metric: w(e) = |v1 - v0|
    Euclidean,
    /// Isotropic conduction velocity map.
    /// Edge weight: w(e) = |e| / harmonic_mean(c[v0], c[v1])
    IsotropicSpeed(Vec<f64>),
    /// Anisotropic fiber-transverse electrophysiological metric.
    ///
    /// At each vertex the diffusion tensor is:
    ///   D(v) = c_f(v)² (f(v)⊗f(v)) + c_t(v)² (I - f(v)⊗f(v))
    ///
    /// The Riemannian metric tensor is g(v) = D(v)^{-1}:
    ///   g(v)(d) = |d|²/c_t(v)² + (d·f(v))² · (1/c_f(v)² − 1/c_t(v)²)
    ///
    /// Edge weight: w(e) = sqrt( (g_v0(d) + g_v1(d)) / 2 )
    /// where d = p(v1) − p(v0).
    FiberTensor {
        /// Per-vertex fiber direction unit vectors (tangent to surface)
        fiber_dirs: Vec<[f64; 3]>,
        /// Per-vertex conduction speed along the fiber direction
        speeds_fiber: Vec<f64>,
        /// Per-vertex conduction speed transverse to the fiber direction
        speeds_transverse: Vec<f64>,
    },
}

#[derive(Debug, Clone)]
pub struct HalfedgeMesh {
    /// One outgoing halfedge per vertex (any, INVALID if isolated)
    pub vert_he: Vec<usize>,
    /// Next halfedge in the same face (CCW)
    pub he_next: Vec<usize>,
    /// Previous halfedge in the same face
    pub he_prev: Vec<usize>,
    /// Twin halfedge (opposite direction on same edge)
    pub he_twin: Vec<usize>,
    /// Source vertex of this halfedge
    pub he_origin: Vec<usize>,
    /// Face to the left of this halfedge (INVALID = boundary)
    pub he_face: Vec<usize>,
    /// Edge index for this halfedge (he and twin share the same edge)
    pub he_edge: Vec<usize>,
    /// One halfedge per face
    pub face_he: Vec<usize>,
    /// 3-D positions
    pub positions: Vec<[f64; 3]>,
    /// Edge lengths indexed by edge index
    pub edge_lengths: Vec<f64>,
    /// Active metric for geodesic computation
    pub metric: Metric,
}

impl HalfedgeMesh {
    /// True when the active metric is plain Euclidean (no speed / fiber map).
    pub fn is_euclidean(&self) -> bool {
        matches!(self.metric, Metric::Euclidean)
    }
    /// Build mesh from vertex positions and triangle faces (0-indexed).
    /// Boundary edges get a twin halfedge with `he_face = INVALID`.
    pub fn from_vertices_faces(verts: &[[f64; 3]], faces: &[[usize; 3]]) -> Self {
        use std::collections::HashMap;

        let n_v = verts.len();
        let n_f = faces.len();

        // --- first pass: allocate interior halfedges 3 per face ---
        let n_interior = n_f * 3;
        let mut he_next = vec![INVALID; n_interior];
        let mut he_prev = vec![INVALID; n_interior];
        let mut he_twin = vec![INVALID; n_interior];
        let mut he_origin = vec![INVALID; n_interior];
        let mut he_face = vec![INVALID; n_interior];
        let mut face_he = vec![INVALID; n_f];

        // directed-edge map: (from, to) -> halfedge index
        let mut dir_map: HashMap<(usize, usize), usize> = HashMap::new();

        for (fi, tri) in faces.iter().enumerate() {
            let base = fi * 3;
            for k in 0..3usize {
                let h = base + k;
                he_origin[h] = tri[k];
                he_face[h] = fi;
                he_next[h] = base + (k + 1) % 3;
                he_prev[h] = base + (k + 2) % 3;
                dir_map.insert((tri[k], tri[(k + 1) % 3]), h);
            }
            face_he[fi] = base;
        }

        // --- link twins for interior edges ---
        for fi in 0..n_f {
            let tri = faces[fi];
            for k in 0..3usize {
                let h = fi * 3 + k;
                let v0 = tri[k];
                let v1 = tri[(k + 1) % 3];
                if let Some(&t) = dir_map.get(&(v1, v0)) {
                    he_twin[h] = t;
                }
            }
        }

        // --- second pass: create boundary halfedges for unmatched edges ---
        let mut bnd_from_h: Vec<usize> = Vec::new(); // interior h that needs a boundary twin
        for h in 0..n_interior {
            if he_twin[h] == INVALID {
                bnd_from_h.push(h);
            }
        }

        let n_bnd = bnd_from_h.len();
        let total_he = n_interior + n_bnd;
        he_next.resize(total_he, INVALID);
        he_prev.resize(total_he, INVALID);
        he_twin.resize(total_he, INVALID);
        he_origin.resize(total_he, INVALID);
        he_face.resize(total_he, INVALID); // boundary face = INVALID

        for (i, &h) in bnd_from_h.iter().enumerate() {
            let bh = n_interior + i;
            // boundary halfedge goes from dest(h) to origin(h)
            he_origin[bh] = he_origin[he_next[h]]; // = dest(h) since he_next[h].origin = v1
            he_face[bh] = INVALID;
            he_twin[h] = bh;
            he_twin[bh] = h;
        }
        // Note: he_next/he_prev for boundary halfedges left INVALID (not needed for our algorithms)

        // --- vert_he: one outgoing halfedge per vertex ---
        let mut vert_he = vec![INVALID; n_v];
        for h in 0..total_he {
            let v = he_origin[h];
            if v < n_v && vert_he[v] == INVALID {
                vert_he[v] = h;
            }
        }

        // --- edge indices: pair each (h, twin) with one edge id ---
        let mut he_edge = vec![INVALID; total_he];
        let mut edge_lengths: Vec<f64> = Vec::new();
        for h in 0..total_he {
            if he_edge[h] == INVALID {
                let t = he_twin[h];
                let eid = edge_lengths.len();
                he_edge[h] = eid;
                if t != INVALID {
                    he_edge[t] = eid;
                }
                let v0 = he_origin[h];
                let v1 = if t != INVALID { he_origin[t] } else { INVALID };
                let len = if v0 < n_v && v1 < n_v {
                    dist3(&verts[v0], &verts[v1])
                } else {
                    0.0
                };
                edge_lengths.push(len);
            }
        }

        HalfedgeMesh {
            vert_he,
            he_next,
            he_prev,
            he_twin,
            he_origin,
            he_face,
            he_edge,
            face_he,
            positions: verts.to_vec(),
            edge_lengths,
            metric: Metric::Euclidean,
        }
    }

    // ---- basic accessors ----

    pub fn n_verts(&self) -> usize { self.positions.len() }
    pub fn n_halfedges(&self) -> usize { self.he_next.len() }
    pub fn n_faces(&self) -> usize { self.face_he.len() }
    pub fn n_edges(&self) -> usize { self.edge_lengths.len() }

    #[inline] pub fn twin(&self, h: usize) -> usize { self.he_twin[h] }
    #[inline] pub fn next(&self, h: usize) -> usize { self.he_next[h] }
    #[inline] pub fn prev(&self, h: usize) -> usize { self.he_prev[h] }
    #[inline] pub fn origin(&self, h: usize) -> usize { self.he_origin[h] }
    #[inline] pub fn dest(&self, h: usize) -> usize { self.he_origin[self.he_twin[h]] }
    #[inline] pub fn face(&self, h: usize) -> usize { self.he_face[h] }
    #[inline] pub fn edge_of(&self, h: usize) -> usize { self.he_edge[h] }
    #[inline] pub fn is_boundary_he(&self, h: usize) -> bool { self.he_face[h] == INVALID }
    #[inline] pub fn is_boundary_edge(&self, h: usize) -> bool {
        self.he_face[h] == INVALID || self.he_face[self.he_twin[h]] == INVALID
    }

    /// Euclidean edge length of edge containing halfedge h
    pub fn edge_len(&self, h: usize) -> f64 {
        self.edge_lengths[self.he_edge[h]]
    }

    /// Metric-weighted edge length according to the active metric mode.
    pub fn metric_edge_len(&self, h: usize) -> f64 {
        let v0 = self.origin(h);
        let v1 = self.dest(h);
        match &self.metric {
            Metric::Euclidean => self.edge_len(h),
            Metric::IsotropicSpeed(speeds) => {
                let base = self.edge_len(h);
                let c = 2.0 / (1.0 / speeds[v0] + 1.0 / speeds[v1]); // harmonic mean speed
                base / c  // travel time = distance / speed
            }
            Metric::FiberTensor { fiber_dirs, speeds_fiber, speeds_transverse } => {
                let p0 = self.positions[v0];
                let p1 = self.positions[v1];
                let d = sub3(&p1, &p0);

                // g_v(d) = |d|²/c_t² + (d·f)²·(1/c_f² − 1/c_t²)
                let quad = |v: usize| {
                    let f  = fiber_dirs[v];
                    let cf = speeds_fiber[v];
                    let ct = speeds_transverse[v];
                    let df = dot3(&d, &f);
                    dot3(&d, &d) / (ct * ct)
                        + df * df * (1.0 / (cf * cf) - 1.0 / (ct * ct))
                };

                // Average metric quadratic form between the two endpoints
                ((quad(v0) + quad(v1)) / 2.0).max(0.0).sqrt()
            }
        }
    }

    /// Set isotropic per-vertex conduction speed map.
    pub fn set_isotropic_speeds(&mut self, speeds: Vec<f64>) {
        self.metric = Metric::IsotropicSpeed(speeds);
    }

    /// Set anisotropic fiber-transverse metric.
    pub fn set_fiber_metric(
        &mut self,
        fiber_dirs: Vec<[f64; 3]>,
        speeds_fiber: Vec<f64>,
        speeds_transverse: Vec<f64>,
    ) {
        self.metric = Metric::FiberTensor { fiber_dirs, speeds_fiber, speeds_transverse };
    }

    /// Reset to Euclidean metric.
    pub fn clear_metric(&mut self) {
        self.metric = Metric::Euclidean;
    }

    /// Outgoing halfedges from vertex v in CCW order
    pub fn outgoing_halfedges(&self, v: usize) -> Vec<usize> {
        let start = self.vert_he[v];
        if start == INVALID { return Vec::new(); }

        let mut result = Vec::new();
        let mut h = start;
        let mut hit_boundary = false;

        // Walk CCW: h → twin(prev(h))
        loop {
            result.push(h);
            let p = self.he_prev[h];
            if p == INVALID { hit_boundary = true; break; }  // h is a boundary halfedge
            let t = self.he_twin[p];
            if t == INVALID { hit_boundary = true; break; }
            if t == start { break; }
            h = t;
        }

        // For boundary vertices: also walk CW (next(twin(h))) from start to
        // collect halfedges on the other arc.
        if hit_boundary {
            let mut h = start;
            loop {
                let t = self.he_twin[h];
                if t == INVALID { break; }
                let n = self.he_next[t];
                if n == INVALID || n == start { break; }
                if result.contains(&n) { break; }
                h = n;
                result.push(h);
            }
        }

        result
    }

    pub fn position(&self, v: usize) -> [f64; 3] { self.positions[v] }

    /// Sum of angles around vertex v (intrinsic angle sum)
    pub fn vertex_angle_sum(&self, v: usize) -> f64 {
        let mut sum = 0.0;
        for h in self.outgoing_halfedges(v) {
            if !self.is_boundary_he(h) {
                sum += self.corner_angle(h);
            }
        }
        sum
    }

    /// Interior angle at the corner of face f at vertex origin(h),
    /// computed from edge lengths using the law of cosines.
    pub fn corner_angle(&self, h: usize) -> f64 {
        // h goes from v0 to v1; prev(h) goes from v2 to v0
        // angle at v0 between edges v0->v1 and v0->v2
        let a = self.edge_len(h);             // v0->v1
        let b = self.edge_len(self.prev(h));   // v2->v0 (same as v0->v2)
        let c = self.edge_len(self.next(h));   // v1->v2 (opposite to v0)
        // law of cosines: c^2 = a^2 + b^2 - 2ab*cos(angle)
        let cos_angle = (a * a + b * b - c * c) / (2.0 * a * b);
        cos_angle.clamp(-1.0, 1.0).acos()
    }

    /// Interior angle at the corner using metric edge lengths.
    /// Same as `corner_angle` but uses `metric_edge_len` so the unfolding
    /// respects the active metric (isotropic speed / fiber tensor).
    pub fn metric_corner_angle(&self, h: usize) -> f64 {
        let a = self.metric_edge_len(h);
        let b = self.metric_edge_len(self.prev(h));
        let c = self.metric_edge_len(self.next(h));
        if a < 1e-30 || b < 1e-30 { return 0.0; }
        let cos_angle = (a * a + b * b - c * c) / (2.0 * a * b);
        cos_angle.clamp(-1.0, 1.0).acos()
    }

    /// 3-D face normal (unnormalized)
    pub fn face_normal(&self, f: usize) -> [f64; 3] {
        let h0 = self.face_he[f];
        let h1 = self.he_next[h0];
        let p0 = self.positions[self.he_origin[h0]];
        let p1 = self.positions[self.he_origin[h1]];
        let p2 = self.positions[self.he_origin[self.he_next[h1]]];
        cross3(&sub3(&p1, &p0), &sub3(&p2, &p0))
    }

    // ────────────────────────────────────────────────────────────────
    // Intrinsic edge flipping
    // ────────────────────────────────────────────────────────────────

    /// Check the intrinsic Delaunay criterion for the edge containing halfedge h.
    /// Uses `corner_angle` (which reads `edge_lengths`). When working on an
    /// intrinsic copy where `edge_lengths` have been set to metric lengths
    /// via `init_metric_edge_lengths()`, this automatically respects the metric.
    pub fn is_intrinsic_delaunay(&self, h: usize) -> bool {
        if self.is_boundary_edge(h) { return true; }
        let t = self.he_twin[h];
        let alpha = self.corner_angle(self.he_next[self.he_next[h]]);
        let beta  = self.corner_angle(self.he_next[self.he_next[t]]);
        alpha + beta <= std::f64::consts::PI + 1e-10
    }

    /// Replace each entry in `edge_lengths` with the metric-weighted edge
    /// length, then set the metric to Euclidean so that `edge_len()` and
    /// `metric_edge_len()` both return the metric length directly.
    /// This prepares the mesh for intrinsic edge flipping: all subsequent
    /// operations (Delaunay checks, fan unfolding, Dijkstra) use the metric
    /// geometry stored as "intrinsic edge lengths".
    pub fn init_metric_edge_lengths(&mut self) {
        let n_he = self.he_next.len();
        // For each edge, find any halfedge and compute metric_edge_len
        let mut visited = vec![false; self.edge_lengths.len()];
        for h in 0..n_he {
            let eid = self.he_edge[h];
            if eid >= self.edge_lengths.len() || visited[eid] { continue; }
            visited[eid] = true;
            self.edge_lengths[eid] = self.metric_edge_len(h);
        }
        // Set metric to Euclidean so metric_edge_len == edge_len
        self.metric = Metric::Euclidean;
    }

    /// Flip the interior edge containing halfedge `h`.
    ///
    /// Before:
    /// ```text
    ///        c                    c
    ///       / \                  /|\
    ///      /   \                / | \
    ///     / f0  \     =>      /  |  \
    ///    a───h──→b           a  f0 f1 b
    ///     \ f1  /             \  |  /
    ///      \   /               \ | /
    ///       \ /                  \|/
    ///        d                    d
    /// ```
    ///
    /// After: the edge becomes c→d.
    /// Edge length is recomputed from the 2-D unfolded quadrilateral.
    ///
    /// Returns `false` if the edge cannot be flipped (boundary, or would
    /// create a degenerate triangle).
    pub fn flip_edge(&mut self, h: usize) -> bool {
        if self.is_boundary_edge(h) { return false; }

        let t = self.he_twin[h];     // twin of h
        let f0 = self.he_face[h];    // left face
        let f1 = self.he_face[t];    // right face

        // Halfedges around the two faces (CCW):
        //   face f0: h → hn → hp   (h: a→b, hn: b→c, hp: c→a)
        //   face f1: t → tn → tp   (t: b→a, tn: a→d, tp: d→b)
        let hn = self.he_next[h];
        let hp = self.he_prev[h];
        let tn = self.he_next[t];
        let tp = self.he_prev[t];

        let a = self.he_origin[h];   // origin of h
        let b = self.he_origin[t];   // origin of t (= dest of h)
        let c = self.he_origin[hp];  // apex of f0 (= dest of hn)
        let d = self.he_origin[tp];  // apex of f1 (= dest of tn)

        // Prevent flipping if it would create a duplicate edge
        // (c and d are already connected)
        for oh in self.outgoing_halfedges(c) {
            if self.dest(oh) == d { return false; }
        }

        // Compute new intrinsic edge length |c−d| by unfolding the quad into 2-D.
        // Uses Euclidean edge lengths for geometric consistency; the metric
        // Delaunay criterion in is_intrinsic_delaunay() drives WHICH edges to flip.
        let lab = self.edge_len(h);
        let lac = self.edge_len(hp);
        let lbc = self.edge_len(hn);
        let lad = self.edge_len(tn);
        let lbd = self.edge_len(tp);

        // Place triangle (a, b, c):
        //   a = (0, 0),  b = (lab, 0)
        //   c from law of cosines at a
        let cos_a_abc = (lab * lab + lac * lac - lbc * lbc) / (2.0 * lab * lac);
        let sin_a_abc = (1.0 - cos_a_abc * cos_a_abc).max(0.0).sqrt();
        let cx = lac * cos_a_abc;
        let cy = lac * sin_a_abc;

        // Place triangle (a, b, d):
        //   d from law of cosines at a, below x-axis
        let cos_a_abd = (lab * lab + lad * lad - lbd * lbd) / (2.0 * lab * lad);
        let sin_a_abd = (1.0 - cos_a_abd * cos_a_abd).max(0.0).sqrt();
        let dx = lad * cos_a_abd;
        let dy = -(lad * sin_a_abd); // below x-axis

        // New edge length: |c - d|
        let new_len = ((cx - dx) * (cx - dx) + (cy - dy) * (cy - dy)).sqrt();
        if new_len < 1e-15 { return false; } // degenerate

        // --- Update halfedge connectivity ---
        // After flip: h becomes c→d, t becomes d→c
        //   face f0: h(c→d) → tp(d→b) → hn(b→c)
        //   face f1: t(d→c) → hp(c→a) → tn(a→d)

        self.he_origin[h] = c;
        self.he_origin[t] = d;

        // face f0: h, tp, hn
        self.he_next[h]  = tp;
        self.he_next[tp] = hn;
        self.he_next[hn] = h;
        self.he_prev[h]  = hn;
        self.he_prev[tp] = h;
        self.he_prev[hn] = tp;
        self.he_face[tp] = f0;

        // face f1: t, hp, tn
        self.he_next[t]  = hp;
        self.he_next[hp] = tn;
        self.he_next[tn] = t;
        self.he_prev[t]  = tn;
        self.he_prev[hp] = t;
        self.he_prev[tn] = hp;
        self.he_face[hp] = f1;

        // face_he: make sure each face points to a valid halfedge
        self.face_he[f0] = h;
        self.face_he[f1] = t;

        // vert_he: vertices a and b might have lost their outgoing halfedge
        // (h was from a, t was from b, but now h is from c and t is from d)
        if self.vert_he[a] == h { self.vert_he[a] = tn; }
        if self.vert_he[b] == t { self.vert_he[b] = tp; }
        // c and d gain new outgoing halfedges
        self.vert_he[c] = h;
        self.vert_he[d] = t;

        // Update edge length (intrinsic)
        let eid = self.he_edge[h];
        self.edge_lengths[eid] = new_len;

        true
    }

    /// Flip all non-Delaunay edges until the triangulation is intrinsically
    /// Delaunay (or `max_flips` is reached). Returns the number of flips performed.
    pub fn flip_to_delaunay(&mut self, max_flips: usize) -> usize {
        let mut total = 0;
        for _ in 0..max_flips {
            let mut flipped_any = false;
            // Iterate over edges (each edge = pair of halfedges)
            let n_he = self.he_next.len();
            for h in 0..n_he {
                if self.he_face[h] == INVALID { continue; } // skip boundary
                let t = self.he_twin[h];
                if h > t { continue; } // process each edge once
                if self.is_intrinsic_delaunay(h) { continue; }
                if self.flip_edge(h) {
                    flipped_any = true;
                    total += 1;
                }
            }
            if !flipped_any { break; }
        }
        total
    }
}

// ---- small vector helpers ----

pub fn dist3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let dz = b[2] - a[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

pub fn sub3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn add3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

pub fn scale3(s: f64, a: &[f64; 3]) -> [f64; 3] {
    [s * a[0], s * a[1], s * a[2]]
}

pub fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn cross3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub fn norm3(a: &[f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

pub fn normalize3(a: &[f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n > 1e-15 { scale3(1.0 / n, a) } else { *a }
}
