/**
 * Author: black_horse2014, chilli
 * Date: 2019-10-29
 * License: Unknown
 * Source: https://codeforces.com/gym/101673/submission/50481926
 * Description: Calculates the area of the union of $n$ polygons (not necessarily
 * convex). The points within each polygon must be given in CCW order.
 * (Epsilon checks may optionally be added to sideOf/sgn, but shouldn't be needed.)
 * T can be double, ll or Fraction. With ll all side tests are exact for
 * $|x|,|y| < 2^{30}$ and only the ratios and the sum are floating point,
 * so prefer it for integer input (in double the tests are only exact
 * for integers up to about $3\cdot 10^7$). With Fraction everything is exact
 * and a Fraction is returned. The area is summed relative to a
 * polygon vertex, so far-away input does not cancel.
 * Time: $O(N^2)$, where $N$ is the total number of points
 * Status: stress-tested, Submitted on ECNA 2017 Problem A
 */
#pragma once

#include "Point.h"
#include "sideOf.h"

typedef Point<double> P;
template<class P> using PolyUnionT = conditional_t<
	is_integral<decltype(P().x)>::value, double, decltype(P().x)>;
template<class P> PolyUnionT<P> rat(P a, P b) {
	typedef PolyUnionT<P> R;
	return sgn(b.x) ? R(a.x) / R(b.x) : R(a.y) / R(b.y);
}
template<class P>
PolyUnionT<P> polyUnion(vector<vector<P>>& poly) {
	typedef decltype(P().x) T; typedef PolyUnionT<P> R;
	R ret = 0;
	P O;
	for (auto& p : poly) if (sz(p)) O = p[0];
	rep(i,0,sz(poly)) rep(v,0,sz(poly[i])) {
		P A = poly[i][v], B = poly[i][(v + 1) % sz(poly[i])];
		vector<pair<R, int>> segs = {{0, 0}, {1, 0}};
		rep(j,0,sz(poly)) if (i != j) {
			rep(u,0,sz(poly[j])) {
				P C = poly[j][u], D = poly[j][(u + 1) % sz(poly[j])];
				int sc = sideOf(A, B, C), sd = sideOf(A, B, D);
				if (sc != sd) {
					T sa = C.cross(D, A), sb = C.cross(D, B);
					if (min(sc, sd) < 0)
						segs.emplace_back(R(sa) / R(sa - sb), sgn(sc - sd));
				} else if (!sc && !sd && j<i && sgn((B-A).dot(D-C))>0){
					segs.emplace_back(rat(C - A, B - A), 1);
					segs.emplace_back(rat(D - A, B - A), -1);
				}
			}
		}
		sort(all(segs));
		for (auto& s : segs) s.first = min(max(s.first, R(0)), R(1));
		R sum = 0;
		int cnt = segs[0].second;
		rep(j,1,sz(segs)) {
			if (!cnt) sum += segs[j].first - segs[j - 1].first;
			cnt += segs[j].second;
		}
		ret += R((A - O).cross(B - O)) * sum;
	}
	return ret / 2;
}
