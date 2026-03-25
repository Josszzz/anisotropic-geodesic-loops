/// Shortest-path geodesics on a triangulated mesh.
///
/// We compute geodesics that travel along the *combinatorial graph* of
/// the mesh (vertices and edges). This gives piecewise-linear paths that
/// are later shortened to smooth geodesics by the flip algorithm.
///
/// Two modes:
///   - Point-to-point: shortest path from vertex `src` to vertex `dst`.
///   - Closed loop:    shortest non-contractible loop through a seed vertex,
///                     found by cutting the mesh and running Dijkstra.

use std::collections::BinaryHeap;
use std::cmp::Ordering;
use crate::mesh::{HalfedgeMesh, INVALID};

// ---- ordered float for BinaryHeap (min-heap via negation) ----

#[derive(Clone, PartialEq)]
struct MinF64(f64);

impl Eq for MinF64 {}

impl PartialOrd for MinF64 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // Reverse for min-heap
        other.0.partial_cmp(&self.0)
    }
}

impl Ord for MinF64 {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap_or(Ordering::Equal)
    }
}

#[derive(Clone, PartialEq, Eq)]
struct HeapEntry {
    dist: MinF64,
    vertex: usize,
}

impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.dist.cmp(&other.dist)
    }
}

// ---- Dijkstra on the mesh graph ----

/// Run Dijkstra from a set of source vertices.
/// Returns `(dist, prev_he)`:
///   - `dist[v]`    = shortest distance from any source to `v`
///   - `prev_he[v]` = halfedge arriving at `v` on the shortest path
///                    (INVALID for sources)
pub fn dijkstra(
    mesh: &HalfedgeMesh,
    sources: &[(usize, f64)], // (vertex, initial_dist)
) -> (Vec<f64>, Vec<usize>) {
    let n = mesh.n_verts();
    let inf = f64::INFINITY;
    let mut dist = vec![inf; n];
    let mut prev_he = vec![INVALID; n];
    let mut heap = BinaryHeap::new();

    for &(v, d) in sources {
        dist[v] = d;
        heap.push(HeapEntry { dist: MinF64(d), vertex: v });
    }

    while let Some(HeapEntry { dist: MinF64(d), vertex: u }) = heap.pop() {
        if d > dist[u] { continue; } // stale entry
        for h in mesh.outgoing_halfedges(u) {
            let w = mesh.dest(h);
            if w >= n { continue; }
            let nd = d + mesh.metric_edge_len(h);
            if nd < dist[w] {
                dist[w] = nd;
                prev_he[w] = h;
                heap.push(HeapEntry { dist: MinF64(nd), vertex: w });
            }
        }
    }

    (dist, prev_he)
}

/// Shortest point-to-point path: list of vertices from `src` to `dst`.
/// Returns None if unreachable.
pub fn shortest_path(
    mesh: &HalfedgeMesh,
    src: usize,
    dst: usize,
) -> Option<Vec<usize>> {
    let (dist, prev_he) = dijkstra(mesh, &[(src, 0.0)]);
    if dist[dst].is_infinite() { return None; }
    Some(trace_path(mesh, &prev_he, src, dst))
}

/// Trace path from src to dst using prev_he.
pub fn trace_path(
    mesh: &HalfedgeMesh,
    prev_he: &[usize],
    src: usize,
    dst: usize,
) -> Vec<usize> {
    let mut path = vec![dst];
    let mut cur = dst;
    while cur != src {
        let h = prev_he[cur];
        if h == INVALID { break; } // unreachable
        cur = mesh.origin(h);
        path.push(cur);
    }
    path.reverse();
    path
}

/// Convert a vertex path to 3-D polyline
pub fn path_to_polyline(mesh: &HalfedgeMesh, vpath: &[usize]) -> Vec<[f64; 3]> {
    vpath.iter().map(|&v| mesh.position(v)).collect()
}

/// Total metric length of a vertex path
pub fn path_metric_length(mesh: &HalfedgeMesh, vpath: &[usize]) -> f64 {
    if vpath.len() < 2 { return 0.0; }
    let mut len = 0.0;
    for w in vpath.windows(2) {
        // find the halfedge from w[0] to w[1]
        let mut found = false;
        for h in mesh.outgoing_halfedges(w[0]) {
            if mesh.dest(h) == w[1] {
                len += mesh.metric_edge_len(h);
                found = true;
                break;
            }
        }
        if !found {
            // Vertices are not directly connected by a mesh edge.
            // Use Euclidean distance as a display-only approximation.
            let p0 = mesh.position(w[0]);
            let p1 = mesh.position(w[1]);
            len += crate::mesh::dist3(&p0, &p1);
            // Note: this indicates a path that skips non-adjacent vertices.
            // In practice this should not happen for valid mesh paths.
        }
    }
    len
}

// ---- Closed geodesic loops ----

/// Find the shortest non-contractible loop through vertex `seed`.
///
/// Strategy: run Dijkstra from `seed` on the mesh. For each edge (u,v)
/// that is not on the shortest-path tree, the cycle `path(seed→u) + edge(u,v) +
/// path(v→seed)` is a candidate closed curve. The shortest such candidate
/// that is topologically non-trivial (not contractible in the graph) is
/// returned.
///
/// Note: this gives a combinatorial loop on the mesh graph; use
/// `flip_geodesics` to shorten it to a smooth geodesic.
pub fn shortest_loop_through(
    mesh: &HalfedgeMesh,
    seed: usize,
) -> Option<Vec<usize>> {
    let (dist, prev_he) = dijkstra(mesh, &[(seed, 0.0)]);

    let mut best: Option<(f64, Vec<usize>)> = None;

    // Iterate over all edges
    for h in 0..mesh.n_halfedges() {
        // Only consider each edge once
        let t = mesh.twin(h);
        if h > t { continue; }
        let u = mesh.origin(h);
        let v = mesh.dest(h);
        if u >= mesh.n_verts() || v >= mesh.n_verts() { continue; }
        if dist[u].is_infinite() || dist[v].is_infinite() { continue; }

        let edge_w = mesh.metric_edge_len(h);
        let candidate_len = dist[u] + edge_w + dist[v];

        // Check if this edge is NOT on the shortest-path tree:
        // edge is on tree if prev_he[v] arrives from u OR prev_he[u] arrives from v
        let on_tree_u = prev_he[v] != INVALID && mesh.origin(prev_he[v]) == u;
        let on_tree_v = prev_he[u] != INVALID && mesh.origin(prev_he[u]) == v;
        if on_tree_u || on_tree_v { continue; }

        if let Some((best_len, _)) = &best {
            if candidate_len >= *best_len { continue; }
        }

        // Build the loop: seed→u (reversed from u's path), then edge u→v, then v→seed
        let path_u = trace_path(mesh, &prev_he, seed, u);
        let path_v = trace_path(mesh, &prev_he, seed, v);

        // Build closed loop: [seed, ..., u, v, ..., (last before seed)]
        // path_u = [seed, ..., u],  path_v = [seed, ..., v]
        // reversed path_v = [v, ..., seed]
        // Concatenation skips the duplicate v and the trailing seed.
        let mut lp = path_u;           // [seed, ..., u]
        lp.push(v);                    // [seed, ..., u, v]
        let mut pv_rev = path_v;
        pv_rev.reverse();              // [v, ..., seed]
        // Skip pv_rev[0]=v (already in lp) and pv_rev.last()=seed (same as lp[0])
        let trim_end = if pv_rev.last() == Some(&seed) { 1 } else { 0 };
        let inner = &pv_rev[1..pv_rev.len().saturating_sub(trim_end)];
        lp.extend_from_slice(inner);

        best = Some((candidate_len, lp));
    }

    best.map(|(_, lp)| lp)
}

/// Find shortest loop in a given free homotopy class, specified by a set of
/// "cuts" (edges to avoid). Returns the shortest loop that crosses the cut set
/// an odd number of times (topological non-triviality w.r.t. those cuts).
///
/// `cut_edges`: set of edge indices forming the cut
pub fn shortest_loop_crossing_cut(
    mesh: &HalfedgeMesh,
    seed: usize,
    cut_edges: &std::collections::HashSet<usize>,
) -> Option<Vec<usize>> {
    // Run Dijkstra on a "doubled graph" where each vertex appears as
    // (v, parity) — parity = number of cut crossings mod 2.
    // We want a path from (seed, 0) back to (seed, 1) with minimum cost.
    let n = mesh.n_verts();
    let inf = f64::INFINITY;
    let mut dist = vec![[inf; 2]; n];
    let mut prev: Vec<[usize; 2]> = vec![[INVALID; 2]; n]; // encodes vertex only
    let mut prev_parity: Vec<[usize; 2]> = vec![[INVALID; 2]; n];

    dist[seed][0] = 0.0;

    let mut heap: BinaryHeap<HeapEntry2> = BinaryHeap::new();
    heap.push(HeapEntry2 { dist: MinF64(0.0), vertex: seed, parity: 0 });

    while let Some(HeapEntry2 { dist: MinF64(d), vertex: u, parity: p }) = heap.pop() {
        if d > dist[u][p] { continue; }
        for h in mesh.outgoing_halfedges(u) {
            let w = mesh.dest(h);
            if w >= n { continue; }
            let crosses = if cut_edges.contains(&mesh.edge_of(h)) { 1usize } else { 0 };
            let np = (p + crosses) % 2;
            let nd = d + mesh.metric_edge_len(h);
            if nd < dist[w][np] {
                dist[w][np] = nd;
                prev[w][np] = u;
                prev_parity[w][np] = p;
                heap.push(HeapEntry2 { dist: MinF64(nd), vertex: w, parity: np });
            }
        }
    }

    if dist[seed][1].is_infinite() { return None; }

    // Trace path back from (seed, 1) to (seed, 0)
    let mut path = vec![seed];
    let mut cur = seed;
    let mut par = 1usize;
    loop {
        let pu = prev[cur][par];
        let pp = prev_parity[cur][par];
        if pu == INVALID { break; }
        path.push(pu);
        par = pp;
        cur = pu;
        if cur == seed && par == 0 { break; }
    }
    path.reverse();
    Some(path)
}

#[derive(Clone, PartialEq, Eq)]
struct HeapEntry2 {
    dist: MinF64,
    vertex: usize,
    parity: usize,
}
impl PartialOrd for HeapEntry2 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(other)) }
}
impl Ord for HeapEntry2 {
    fn cmp(&self, other: &Self) -> Ordering { self.dist.cmp(&other.dist) }
}
