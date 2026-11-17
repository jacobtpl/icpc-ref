/**
 * Author: Devin
 * Date: 2026-11-16
 * License: CC0
 * Source: randomized incremental hull with conflict lists
 *  (de Berg et al., Computational Geometry, chapter 11)
 * Description: Faces of the 3D convex hull as index triples into A,
 *  counter-clockwise seen from outside. Any Point3D<T>.
 *  d sees face abc if (b-a) x (c-a) . (d-a) > eps.
 *  Integer T: exact with eps = 0 while 6D^3 fits in T, D = largest
 *  coordinate difference (ll: 1e6, __int128: 2e12). Duplicate, collinear
 *  and coplanar points are fine; coplanar facets come out triangulated.
 *  Floating T: eps = 0 in general position; eps = 1e-9 * D^3 (double) if
 *  the input is degenerate up to rounding (e.g. a grid times 0.1).
 *  Points a little off a plane (eps/10 to 1e4*eps) can end up outside
 *  the result for any eps: scale to integers instead.
 *  The result is always a closed surface. Asserts if all points are
 *  coplanar. On a hackable judge seed the shuffle from the clock.
 * Time: expected O(n log n)
 * Memory: expected O(n log n)
 * Usage: hull3dFast(A) or hull3dFast(A, 1e-9 * D*D*D)
 * Status: stress-tested against brute force and 3dHull.h
 */
#pragma once

#include "Point3D.h"

template<class P3>
vector<array<int,3>> hull3dFast(const vector<P3>& A,
		decltype(P3().x) eps = 0) {
	int n = sz(A);
	vi o(n); iota(all(o), 0); shuffle(all(o), mt19937(5));
	auto vol = [&](int a, int b, int c, int d) { // >eps: d sees abc
		return (A[b]-A[a]).cross(A[c]-A[a]).dot(A[d]-A[a]); };
	auto ab = [](auto x) { return x < 0 ? -x : x; };
	// start from a fat tetrahedron: farthest point, then farthest
	// from that line, then farthest from that plane
	auto far = [&](int k, auto f) {
		rep(i,k+1,n) if (f(o[k]) < f(o[i])) swap(o[k], o[i]); };
	far(1, [&](int i) { return (A[i]-A[o[0]]).dist2(); });
	far(2, [&](int i) {
		P3 c = (A[o[1]]-A[o[0]]).cross(A[i]-A[o[0]]);
		return max({ab(c.x), ab(c.y), ab(c.z)}); });
	far(3, [&](int i) { return ab(vol(o[0], o[1], o[2], i)); });
	assert(n >= 4 && ab(vol(o[0], o[1], o[2], o[3])) > eps);
	if (vol(o[0], o[1], o[2], o[3]) > 0) swap(o[1], o[2]);

	vector<array<int,3>> fs, nb, res; // faces ever made; face across edge j,j+1
	vector<bool> dead;
	vector<vi> pts, vis(n); // points seeing face, faces seen by point
	vi seen(n, -1), at(n), in(n, -1);
	auto add = [&](int a, int b, int c, const vi& cand) {
		int id = sz(fs);
		fs.push_back({a, b, c}); dead.push_back(0);
		vi ps;
		for (int p : cand) if (seen[p] != id) {
			seen[p] = id;
			if (p != c && vol(a, b, c, p) > eps)
				ps.push_back(p), vis[p].push_back(id);
		}
		pts.push_back(ps);
	};
	vi rest(o.begin() + 4, o.end());
	add(o[0], o[1], o[2], rest); add(o[0], o[3], o[1], rest);
	add(o[1], o[3], o[2], rest); add(o[2], o[3], o[0], rest);
	nb = {{1, 2, 3}, {3, 2, 0}, {1, 3, 0}, {2, 1, 0}};
	for (int i : rest) {
		vi rem; // faces seen from i: replace them by a cone from i.
		// Grow them from one face, only in steps that keep them a disc
		// (two shared edges, or one and a new corner), so the horizon
		// is one simple cycle whatever the tests say
		auto kill = [&](int f) {
			dead[f] = 1, rem.push_back(f);
			for (int v : fs[f]) in[v] = i;
		};
		for (int f : vis[i]) if (!dead[f]) { kill(f); break; }
		// second pass: a face i clearly sees is stuck, so also take
		// faces i is within eps of, to fill the pocket that blocks it
		rep(pass,0,2) {
			bool stuck = 0;
			rep(q,0,sz(rem)) for (int g : nb[rem[q]]) if (!dead[g]) {
				auto v = vol(fs[g][0], fs[g][1], fs[g][2], i);
				if (!(v > (pass ? -eps : eps))) continue;
				int c = 0, w = 0;
				rep(k,0,3) if (dead[nb[g][k]]) c++, w = fs[g][(k+2)%3];
				if (c == 2 || (c == 1 && in[w] != i)) kill(g);
				else if (v > eps) stuck = 1;
			}
			if (!stuck) break;
		}
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
