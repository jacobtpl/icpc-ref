/**
 * Author: Devin
 * Date: 2026-11-16
 * License: CC0
 * Source: randomized incremental hull with conflict lists
 *  (de Berg et al., Computational Geometry, chapter 11)
 * Description: Computes all faces of the 3-dimensional hull of a point set
 *  as index triples into A, counter-clockwise seen from outside.
 *  Works on any Point3D<T>. With exact arithmetic (integer T) coplanar,
 *  collinear and duplicate points are fine: coplanar facets come out
 *  triangulated and points that are not corners of the hull may or may
 *  not be used. Only fails (assert) if all points are coplanar.
 *  Every test is vol() on coordinate differences, so only the extent D
 *  (largest coordinate difference) matters: exact while 6D^3 fits in T
 *  (long long: D <= 1e6, __int128: D <= 2e12, double with
 *  integer coordinates: D <= 1.1e5). Non-integer doubles need
 *  general position.
 *  For a hackable judge seed the shuffle from the clock.
 * Time: expected O(n log n)
 * Memory: expected O(n log n)
 * Status: stress-tested against brute force and 3dHull.h
 */
#pragma once

#include "Point3D.h"

template<class P3>
vector<array<int,3>> hull3dFast(const vector<P3>& A) {
	int n = sz(A);
	vi o(n); iota(all(o), 0); shuffle(all(o), mt19937(5));
	auto vol = [&](int a, int b, int c, int d) { // >0: d sees abc
		return (A[b]-A[a]).cross(A[c]-A[a]).dot(A[d]-A[a]); };
	// move a non-degenerate tetrahedron to the front
	rep(i,1,n) if (A[o[i]] != A[o[0]]) { swap(o[1], o[i]); break; }
	rep(i,2,n) if ((A[o[1]]-A[o[0]]).cross(A[o[i]]-A[o[0]])
		!= P3()) { swap(o[2], o[i]); break; }
	rep(i,3,n) if (vol(o[0], o[1], o[2], o[i]) != 0) {
		swap(o[3], o[i]); break; }
	assert(n >= 4 && vol(o[0], o[1], o[2], o[3]) != 0);
	if (vol(o[0], o[1], o[2], o[3]) > 0) swap(o[1], o[2]);

	vector<array<int,3>> fs, nb, res; // faces ever made; face across edge j,j+1
	vector<bool> dead;
	vector<vi> pts, vis(n); // points seeing face, faces seen by point
	vi seen(n, -1), at(n);
	auto add = [&](int a, int b, int c, const vi& cand) {
		int id = sz(fs);
		fs.push_back({a, b, c}); dead.push_back(0);
		vi ps;
		for (int p : cand) if (seen[p] != id) {
			seen[p] = id;
			if (p != c && vol(a, b, c, p) > 0)
				ps.push_back(p), vis[p].push_back(id);
		}
		pts.push_back(ps);
	};
	vi rest(o.begin() + 4, o.end());
	add(o[0], o[1], o[2], rest); add(o[0], o[3], o[1], rest);
	add(o[1], o[3], o[2], rest); add(o[2], o[3], o[0], rest);
	nb = {{1, 2, 3}, {3, 2, 0}, {1, 3, 0}, {2, 1, 0}};
	for (int i : rest) {
		vi rem; // faces seen from i: replace them by a cone from i
		for (int f : vis[i]) if (!dead[f]) dead[f] = 1, rem.push_back(f);
		int st = sz(fs);
		for (int f : rem) rep(j,0,3) {
			int a = fs[f][j], b = fs[f][(j+1)%3], g = nb[f][j];
			if (dead[g]) continue; // else a->b is on the horizon
			vi cand = pts[f];
			cand.insert(cand.end(), all(pts[g]));
			at[a] = sz(fs);
			nb.push_back({g, 0, 0});
			rep(k,0,3) if (fs[g][k] == b) nb[g][k] = sz(fs);
			add(a, b, i, cand);
		}
		rep(f,st,sz(fs)) { // glue the cone: face b,i meets at[b]
			int h = at[fs[f][1]];
			nb[f][1] = h, nb[h][2] = f;
		}
		for (int f : rem) vi().swap(pts[f]);
	}
	rep(f,0,sz(fs)) if (!dead[f]) res.push_back(fs[f]);
	return res;
}
