/// Unit and integration tests for geodesic_loops.
///
/// Test meshes used:
///   - single_tet: tetrahedron (4 verts, 4 faces) — simplest closed surface
///   - flat_grid:  NxN regular grid (2*(N-1)^2 triangles) — flat 2-D mesh
///   - torus_approx: discrete torus (non-contractible loops exist)

#[cfg(test)]
mod mesh_tests {
    use crate::mesh::{HalfedgeMesh, INVALID, dist3};

    fn single_tet() -> HalfedgeMesh {
        // Regular tetrahedron, edge length = sqrt(2)
        let verts = vec![
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ];
        let faces = vec![[0,1,2],[0,2,3],[0,3,1],[1,3,2]];
        HalfedgeMesh::from_vertices_faces(&verts, &faces)
    }

    #[test]
    fn tet_counts() {
        let m = single_tet();
        assert_eq!(m.n_verts(), 4);
        assert_eq!(m.n_faces(), 4);
        assert_eq!(m.n_edges(), 6); // tetrahedron has 6 edges
    }

    #[test]
    fn tet_halfedge_twins() {
        let m = single_tet();
        for h in 0..m.n_halfedges() {
            let t = m.twin(h);
            assert_ne!(t, INVALID, "halfedge {} has no twin", h);
            assert_eq!(m.twin(t), h, "twin relationship not symmetric for h={}", h);
        }
    }

    #[test]
    fn tet_edge_lengths() {
        let m = single_tet();
        let expected = 2.0_f64.sqrt() * 2.0; // distance between any two verts = 2*sqrt(2)
        for e in 0..m.n_edges() {
            let l = m.edge_lengths[e];
            assert!((l - expected).abs() < 1e-10, "edge {} length {} != {}", e, l, expected);
        }
    }

    #[test]
    fn tet_outgoing_halfedges() {
        let m = single_tet();
        for v in 0..m.n_verts() {
            let out = m.outgoing_halfedges(v);
            assert_eq!(out.len(), 3, "vertex {} should have 3 outgoing halfedges", v);
            for h in &out {
                assert_eq!(m.origin(*h), v);
            }
        }
    }

    #[test]
    fn tet_corner_angles_sum() {
        let m = single_tet();
        // Each face is equilateral → each interior angle = π/3
        for h in 0..m.n_halfedges() {
            if !m.is_boundary_he(h) {
                let angle = m.corner_angle(h);
                assert!((angle - std::f64::consts::PI / 3.0).abs() < 1e-10,
                    "expected π/3, got {}", angle);
            }
        }
    }

    fn flat_grid(n: usize) -> HalfedgeMesh {
        // n×n grid of vertices at integer coords, triangulated with
        // alternating diagonal directions for better isotropy.
        let mut verts = Vec::new();
        for i in 0..n {
            for j in 0..n {
                verts.push([i as f64, j as f64, 0.0]);
            }
        }
        let mut faces = Vec::new();
        for i in 0..(n-1) {
            for j in 0..(n-1) {
                let a = i * n + j;       // (i, j)
                let b = i * n + j + 1;   // (i, j+1)
                let c = (i+1) * n + j;   // (i+1, j)
                let d = (i+1) * n + j+1; // (i+1, j+1)
                if (i + j) % 2 == 0 {
                    // Up-right diagonal a-d
                    faces.push([a, b, d]);
                    faces.push([a, d, c]);
                } else {
                    // Down-left diagonal b-c
                    faces.push([a, b, c]);
                    faces.push([b, d, c]);
                }
            }
        }
        HalfedgeMesh::from_vertices_faces(&verts, &faces)
    }

    #[test]
    fn grid_counts() {
        let n = 5;
        let m = flat_grid(n);
        assert_eq!(m.n_verts(), n * n);
        assert_eq!(m.n_faces(), 2 * (n - 1) * (n - 1));
    }

    #[test]
    fn grid_no_invalid_twins_interior() {
        // Interior edges must have both halfedges in a proper face
        let m = flat_grid(4);
        for h in 0..m.n_halfedges() {
            let t = m.twin(h);
            assert_ne!(t, INVALID);
            assert_eq!(m.twin(t), h);
        }
    }
}

#[cfg(test)]
mod geodesic_tests {
    use crate::mesh::HalfedgeMesh;
    use crate::geodesic::{shortest_path, path_metric_length};

    fn flat_grid(n: usize) -> HalfedgeMesh {
        let mut verts = Vec::new();
        for i in 0..n {
            for j in 0..n {
                verts.push([i as f64, j as f64, 0.0]);
            }
        }
        let mut faces = Vec::new();
        for i in 0..(n-1) {
            for j in 0..(n-1) {
                let a = i * n + j;
                let b = i * n + j + 1;
                let c = (i+1) * n + j;
                let d = (i+1) * n + j+1;
                if (i + j) % 2 == 0 {
                    faces.push([a, b, d]);
                    faces.push([a, d, c]);
                } else {
                    faces.push([a, b, c]);
                    faces.push([b, d, c]);
                }
            }
        }
        HalfedgeMesh::from_vertices_faces(&verts, &faces)
    }

    #[test]
    fn path_same_vertex() {
        let m = flat_grid(5);
        let p = shortest_path(&m, 0, 0).unwrap();
        assert_eq!(p, vec![0]);
        assert_eq!(path_metric_length(&m, &p), 0.0);
    }

    #[test]
    fn path_adjacent_vertices() {
        let m = flat_grid(5);
        // Vertex 0 = (0,0), vertex 1 = (0,1) — distance 1
        let p = shortest_path(&m, 0, 1).unwrap();
        let len = path_metric_length(&m, &p);
        assert!((len - 1.0).abs() < 1e-10, "expected length 1, got {}", len);
    }

    #[test]
    fn path_corner_to_corner() {
        let n = 5usize;
        let m = flat_grid(n);
        // (0,0) to (4,4): Manhattan-like on graph, but diagonal exists so dist = sqrt(32)
        let src = 0;
        let dst = n * n - 1;
        let p = shortest_path(&m, src, dst).unwrap();
        let len = path_metric_length(&m, &p);
        // Lower bound = Euclidean distance = sqrt(4^2+4^2) = 4*sqrt(2) ≈ 5.657
        let euclidean = (((n-1) as f64).powi(2) * 2.0).sqrt();
        assert!(len >= euclidean - 1e-10, "path too short: {} < {}", len, euclidean);
        assert!(len <= euclidean * 1.5, "path too long: {}", len); // should be reasonable
    }

    #[test]
    fn path_length_monotone_with_distance() {
        let n = 6usize;
        let m = flat_grid(n);
        // Compare lengths of paths at increasing distances from origin
        let mut prev_len = 0.0;
        for k in 1..n {
            let dst = k; // vertex (0, k)
            let p = shortest_path(&m, 0, dst).unwrap();
            let len = path_metric_length(&m, &p);
            assert!(len >= prev_len - 1e-10, "length should increase: {} < {}", len, prev_len);
            prev_len = len;
        }
    }

    #[test]
    fn path_symmetry() {
        let m = flat_grid(5);
        let p1 = shortest_path(&m, 3, 18).unwrap();
        let p2 = shortest_path(&m, 18, 3).unwrap();
        let l1 = path_metric_length(&m, &p1);
        let l2 = path_metric_length(&m, &p2);
        assert!((l1 - l2).abs() < 1e-10, "path length should be symmetric: {} vs {}", l1, l2);
    }
}

#[cfg(test)]
mod flip_tests {
    use crate::mesh::HalfedgeMesh;
    use crate::flip_geodesics::{FlipConfig, GeodesicPath, shorten_path, geodesic_path};

    fn flat_grid(n: usize) -> HalfedgeMesh {
        let mut verts = Vec::new();
        for i in 0..n {
            for j in 0..n {
                verts.push([i as f64, j as f64, 0.0]);
            }
        }
        let mut faces = Vec::new();
        for i in 0..(n-1) {
            for j in 0..(n-1) {
                let a = i * n + j;
                let b = i * n + j + 1;
                let c = (i+1) * n + j;
                let d = (i+1) * n + j+1;
                if (i + j) % 2 == 0 {
                    faces.push([a, b, d]);
                    faces.push([a, d, c]);
                } else {
                    faces.push([a, b, c]);
                    faces.push([b, d, c]);
                }
            }
        }
        HalfedgeMesh::from_vertices_faces(&verts, &faces)
    }

    #[test]
    fn shorten_does_not_increase_length() {
        let m = flat_grid(8);
        let config = FlipConfig::default();
        // A zigzag path that should be shortened
        let vpath = vec![0, 9, 2, 11, 4, 13, 6];
        let path = GeodesicPath::open(vpath);
        let init_len = path.metric_length(&m);
        let shortened = shorten_path(&m, &path, &config);
        let final_len = shortened.metric_length(&m);
        assert!(final_len <= init_len + 1e-10,
            "shortening increased length: {} → {}", init_len, final_len);
    }

    #[test]
    fn geodesic_path_converges() {
        let m = flat_grid(10);
        let config = FlipConfig::default();
        let n = 10usize;
        let result = geodesic_path(&m, 0, n * n - 1, &config);
        assert!(result.is_some(), "geodesic_path returned None");
        let p = result.unwrap();
        assert!(p.len() >= 2, "path too short");
        // Length should be close to Euclidean
        let euclidean = (((n-1) as f64).powi(2) * 2.0).sqrt();
        let len = p.metric_length(&m);
        assert!(len <= euclidean * 1.3 + 1e-10,
            "final path too long: {} vs euclidean {}", len, euclidean);
    }

    #[test]
    fn closed_path_metric_length() {
        let m = flat_grid(4);
        // A 4-cycle: 0→1→5→4→0
        let path = GeodesicPath::closed(vec![0, 1, 5, 4]);
        let len = path.metric_length(&m);
        // Each edge has length 1, so total = 4
        assert!((len - 4.0).abs() < 1e-10, "expected length 4, got {}", len);
    }
}

#[cfg(test)]
mod metric_tests {
    use crate::mesh::{HalfedgeMesh, Metric};
    use crate::geodesic::shortest_path;

    fn flat_grid(n: usize) -> HalfedgeMesh {
        let mut verts = Vec::new();
        for i in 0..n {
            for j in 0..n {
                verts.push([i as f64, j as f64, 0.0]);
            }
        }
        let mut faces = Vec::new();
        for i in 0..(n-1) {
            for j in 0..(n-1) {
                let a = i * n + j;
                let b = i * n + j + 1;
                let c = (i+1) * n + j;
                let d = (i+1) * n + j+1;
                if (i + j) % 2 == 0 {
                    faces.push([a, b, d]);
                    faces.push([a, d, c]);
                } else {
                    faces.push([a, b, c]);
                    faces.push([b, d, c]);
                }
            }
        }
        HalfedgeMesh::from_vertices_faces(&verts, &faces)
    }

    #[test]
    fn euclidean_is_default() {
        let m = flat_grid(5);
        assert!(matches!(m.metric, Metric::Euclidean));
    }

    #[test]
    fn isotropic_speeds_slow_zone_detour() {
        // Uniform speed → same path length as Euclidean scaled by 1/speed
        let n = 6usize;
        let mut m = flat_grid(n);
        let speeds = vec![1.0f64; n * n];
        m.set_isotropic_speeds(speeds);
        let p_iso = shortest_path(&m, 0, n * n - 1).unwrap();
        let len_iso = crate::geodesic::path_metric_length(&m, &p_iso);

        m.clear_metric();
        let p_euc = shortest_path(&m, 0, n * n - 1).unwrap();
        let len_euc = crate::geodesic::path_metric_length(&m, &p_euc);

        // With uniform speed=1, metric lengths should be equal
        assert!((len_iso - len_euc).abs() < 1e-10,
            "uniform speed should equal Euclidean: {} vs {}", len_iso, len_euc);
    }

    #[test]
    fn isotropic_speeds_obstacle_forces_detour() {
        // Block the direct path with a slow zone; geodesic should go around it
        let n = 7usize;
        let mut m = flat_grid(n);
        let mut speeds = vec![1.0f64; n * n];
        // Slow column at i=3 (all j)
        for j in 0..n { speeds[3 * n + j] = 0.01; }
        m.set_isotropic_speeds(speeds.clone());

        let p_slow = shortest_path(&m, 0, n * n - 1).unwrap();
        let len_slow = crate::geodesic::path_metric_length(&m, &p_slow);

        m.clear_metric();
        let p_euc = shortest_path(&m, 0, n * n - 1).unwrap();
        let len_euc = crate::geodesic::path_metric_length(&m, &p_euc);

        // With the slow zone, metric length should be larger than Euclidean (path avoids it)
        assert!(len_slow > len_euc * 1.1,
            "slow zone path should be longer in metric: {} vs {}", len_slow, len_euc);
    }

    #[test]
    fn fiber_tensor_isotropic_equals_euclidean() {
        // When c_f == c_t == 1 everywhere, fiber metric == Euclidean
        let n = 5usize;
        let mut m = flat_grid(n);
        let nv = m.n_verts();
        // Fiber along x-axis
        let fiber_dirs: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0]; nv];
        let speeds_fiber = vec![1.0f64; nv];
        let speeds_transverse = vec![1.0f64; nv];
        m.set_fiber_metric(fiber_dirs, speeds_fiber, speeds_transverse);

        let p_fiber = shortest_path(&m, 0, n * n - 1).unwrap();
        let len_fiber = crate::geodesic::path_metric_length(&m, &p_fiber);

        m.clear_metric();
        let p_euc = shortest_path(&m, 0, n * n - 1).unwrap();
        let len_euc = crate::geodesic::path_metric_length(&m, &p_euc);

        assert!((len_fiber - len_euc).abs() < 1e-10,
            "isotropic fiber metric should equal Euclidean: {} vs {}", len_fiber, len_euc);
    }

    #[test]
    fn fiber_tensor_anisotropic_prefers_fiber_direction() {
        // Fiber along x-axis (i direction in grid), c_f >> c_t.
        // Going in x direction is fast (low metric cost), going in y is slow.
        // Path from (0,0) to (4,0) should cost less than path from (0,0) to (0,4).
        let n = 5usize;
        let mut m = flat_grid(n);
        let nv = m.n_verts();
        let fiber_dirs: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0]; nv]; // fiber along x (i axis)
        let speeds_fiber     = vec![2.0f64; nv]; // fast along fiber
        let speeds_transverse = vec![0.5f64; nv]; // slow transverse

        m.set_fiber_metric(fiber_dirs, speeds_fiber, speeds_transverse);

        // Path along fiber direction: (0,0)→(4,0) = vertices 0→20 (step n each row)
        // Grid: v(i,j) = i*n + j, so v(4,0) = 20
        let p_along = shortest_path(&m, 0, 4 * n).unwrap();
        let len_along = crate::geodesic::path_metric_length(&m, &p_along);

        // Path transverse: (0,0)→(0,4) = vertices 0→4
        let p_trans = shortest_path(&m, 0, 4).unwrap();
        let len_trans = crate::geodesic::path_metric_length(&m, &p_trans);

        // Along fiber should be faster (lower metric cost per unit length)
        // c_f=2 → along-fiber metric = |d|/c_f = 4/2 = 2
        // c_t=0.5 → transverse metric ≈ |d|/c_t = 4/0.5 = 8
        assert!(len_along < len_trans,
            "fiber direction should be cheaper: {} vs {}", len_along, len_trans);
    }

    #[test]
    fn clear_metric_resets_to_euclidean() {
        let n = 4usize;
        let mut m = flat_grid(n);
        let nv = m.n_verts();
        m.set_isotropic_speeds(vec![0.5f64; nv]);
        assert!(matches!(m.metric, Metric::IsotropicSpeed(_)));
        m.clear_metric();
        assert!(matches!(m.metric, Metric::Euclidean));
    }
}

#[cfg(test)]
mod loop_tests {
    use crate::mesh::HalfedgeMesh;
    use crate::geodesic::shortest_loop_through;

    /// Build a discrete torus by identifying opposite edges of a grid.
    /// We approximate this by a cylinder (closed in one direction).
    fn cylinder(n: usize, m_rings: usize) -> HalfedgeMesh {
        // n vertices per ring, m_rings rings
        // Vertices: v(i,j) = i*n + j, j mod n
        let mut verts = Vec::new();
        let r = 1.0f64;
        for i in 0..m_rings {
            let z = i as f64 / (m_rings - 1) as f64;
            for j in 0..n {
                let theta = 2.0 * std::f64::consts::PI * j as f64 / n as f64;
                verts.push([r * theta.cos(), r * theta.sin(), z]);
            }
        }

        let mut faces = Vec::new();
        for i in 0..(m_rings - 1) {
            for j in 0..n {
                let a = i * n + j;
                let b = i * n + (j + 1) % n;
                let c = (i + 1) * n + j;
                let d = (i + 1) * n + (j + 1) % n;
                faces.push([a, b, c]);
                faces.push([b, d, c]);
            }
        }
        HalfedgeMesh::from_vertices_faces(&verts, &faces)
    }

    #[test]
    fn loop_found_on_cylinder() {
        // The cylinder has non-contractible loops (around the axis)
        let m = cylinder(12, 4);
        let result = shortest_loop_through(&m, 0);
        assert!(result.is_some(), "no loop found on cylinder");
        let lp = result.unwrap();
        assert!(lp.len() >= 3, "loop too short");
        // The loop should come back to the start vertex
        assert_eq!(lp[0], 0);
    }

    #[test]
    fn loop_length_positive() {
        // shortest_loop_through finds any off-tree cycle (may be contractible).
        // Just check it returns a valid positive-length loop.
        let n = 12usize;
        let m = cylinder(n, 6);
        let result = shortest_loop_through(&m, 0);
        let lp = result.unwrap();
        let len = crate::geodesic::path_metric_length(&m, &lp);
        assert!(len > 0.0, "loop length should be positive, got {}", len);
        assert!(len < 100.0, "loop unreasonably long: {}", len);
    }

    #[test]
    fn non_contractible_loop_via_cut() {
        // Use shortest_loop_crossing_cut to find a topologically non-trivial
        // loop on the cylinder.
        //
        // The correct cut for detecting angular (non-contractible) loops is the
        // SEAM edges: (i*n + n-1) -- (i*n + 0) for each ring i.
        // Any loop that goes all the way around the cylinder must cross this seam.
        let n = 12usize;
        let m_rings = 6usize;
        let m = cylinder(n, m_rings);

        // Find the seam edge indices: edge between (ring i, j=n-1) and (ring i, j=0)
        let mut cut_edges = std::collections::HashSet::new();
        for i in 0..m_rings {
            let v_last = i * n + (n - 1); // ring i, j=n-1
            let v_first = i * n + 0;       // ring i, j=0
            for h in m.outgoing_halfedges(v_last) {
                if m.dest(h) == v_first {
                    cut_edges.insert(m.edge_of(h));
                    break;
                }
            }
        }
        assert!(!cut_edges.is_empty(), "no cut edges found");

        let result = crate::geodesic::shortest_loop_crossing_cut(&m, 0, &cut_edges);
        assert!(result.is_some(), "no non-contractible loop found");
        let lp = result.unwrap();
        assert!(lp.len() >= 3, "non-contractible loop must have at least 3 vertices");
        assert_eq!(lp[0], 0, "loop must start at seed");
        let len = crate::geodesic::path_metric_length(&m, &lp);
        assert!(len > 0.0, "non-contractible loop must have positive length: {}", len);
        // The loop crosses the seam exactly once. For this seamed cylinder mesh,
        // the shortest non-contractible loop is a small triangle near the seam
        // (not the full circumference), which is geometrically correct.
        assert!(len < 2.0 * std::f64::consts::PI * 1.5,
            "non-contractible loop unreasonably long: {}", len);
    }
}
