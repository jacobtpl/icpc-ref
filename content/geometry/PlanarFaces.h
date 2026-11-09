/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: own work; standard half-edge face tracing (de Berg et al., Computational Geometry, ch. 2)
 * Description: Faces of a planar straight-line graph. p are distinct points,
 *  e are distinct edges (0-indexed endpoints, no loops) that only meet at common endpoints.
 *  Half-edge $2i$ is e[i].first $\to$ e[i].second, $2i+1$ its reverse; fid[h] is the face left of h.
 *  f[k] lists the vertices of boundary walk $k$ in order (a vertex repeats at bridges/cut vertices),
 *  ar[k] is twice its signed area. Bounded faces are CCW with ar $> 0$; every component
 *  with an edge gives exactly one outer walk, CW with ar $\le 0$ ($0$ iff it is a tree).
 *  Isolated vertices give nothing. Components are not nested into each other: a face with
 *  holes appears as its CCW walk plus the outer walk of every component inside.
 *  Exact for $|x|, |y| \le 10^9$.
 * Time: O((n + m) \log m)
 * Usage: PlanarFaces F(p, e); rep(k,0,sz(F.f)) if (F.ar[k] > 0) ...
 * Status: stress-tested against a brute-force tracer and exact winding numbers
 */
#pragma once

#include "Point.h"

typedef Point<ll> P;
struct PlanarFaces {
	vi fid; vector<vi> f; vector<ll> ar;
	PlanarFaces(const vector<P>& p, const vector<pii>& e) {
		int m = 2 * sz(e);
		auto to = [&](int h) {
			return h & 1 ? e[h / 2].first : e[h / 2].second; };
		auto dir = [&](int h) { return p[to(h)] - p[to(h ^ 1)]; };
		auto half = [](P d) {
			return d.y < 0 || (d.y == 0 && d.x < 0); };
		vector<vi> g(sz(p));
		rep(h,0,m) g[to(h ^ 1)].push_back(h);
		vi nx(m);
		for (vi& v : g) {
			sort(all(v), [&](int a, int b) {
				P c = dir(a), d = dir(b);
				return half(c) != half(d) ? half(c) < half(d)
				                          : c.cross(d) > 0;
			});
			rep(j,0,sz(v)) nx[v[j] ^ 1] = v[j ? j - 1 : sz(v) - 1];
		}
		fid.assign(m, -1);
		rep(i,0,m) if (fid[i] < 0) {
			unsigned long long a = 0; // wraps, final value fits
			f.emplace_back();
			for (int h = i; fid[h] < 0; h = nx[h]) {
				fid[h] = sz(f) - 1;
				f.back().push_back(to(h ^ 1));
				a += p[to(h ^ 1)].cross(p[to(h)]);
			}
			ar.push_back((ll)a);
		}
	}
};
