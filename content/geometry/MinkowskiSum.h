/**
 * Author: Benq
 * Description: Minkowski sum of two convex polygons given in CCW order.
 * Either input may also be a single point, a segment (two points) or
 * have collinear/repeated boundary vertices; an empty input gives \{\}.
 * Returns the strictly convex hull of the sum in CCW order starting
 * from its lexicographically smallest vertex (1 or 2 points if degenerate).
 * P is Point<T>; with T = ll it is exact for $|x|,|y| < 2^{30}$. Output
 * coordinates are up to twice the input, so halve that bound if the
 * result is fed to another cross product.
 * example: diameter of hull A is max in A+(-A)
 * Time: O(N)
 * Source: https://cp-algorithms.com/geometry/minkowski.html#implementation
 * Status: stress-tested against the hull of all pairwise sums; tested on
   * MITPC
   * https://open.kattis.com/problems/roberthood
 */
#pragma once

#include "Point.h"

template<class P>
vector<P> minkowski_sum(vector<P> a, vector<P> b) {
	if (a.empty() || b.empty()) return {};
	auto edges = [](vector<P>& v) {
		rotate(begin(v), min_element(all(v)), end(v));
		vector<P> e;
		rep(i,0,sz(v)) {
			P d = v[(i + 1) % sz(v)] - v[i];
			if (!(d == P())) e.push_back(d);
		}
		return e; // sorted by angle in (-pi/2, 3pi/2]
	};
	auto up = [](P d) { return d.x > 0 || (d.x == 0 && d.y > 0); };
	auto cmp = [&](P u, P v) {
		return up(u) != up(v) ? up(u) : u.cross(v) > 0; };
	vector<P> ea = edges(a), eb = edges(b), e(sz(ea) + sz(eb));
	merge(all(ea), all(eb), begin(e), cmp);
	vector<P> res = {a[0] + b[0]};
	rep(i,0,sz(e)) {
		P p = res.back() + e[i];
		if (i && !cmp(e[i-1], e[i])) res.back() = p; // parallel
		else res.push_back(p);
	}
	if (sz(res) > 1) res.pop_back(); // back at the start
	return res;
}
