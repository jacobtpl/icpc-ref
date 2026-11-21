/**
 * Author: black_horse2014, chilli
 * Date: 2019-10-29
 * License: Unknown
 * Source: https://codeforces.com/gym/101673/submission/50481926
 * Description: Calculates the area of the union of $n$ polygons (not necessarily
 * convex). The points within each polygon must be given in CCW order.
 * Returns double for int and ll input, else T. With integer T all side
 * tests are exact (ll: |x|,|y| < 2^30, int: < 2^14), so prefer it for
 * integer input (in double the tests are only exact up to about 3e7).
 * Eps: side tests are cross products through sgn, so geoEps applies;
 * it shouldn't be needed. The area is summed relative to a
 * polygon vertex O, so far-away input does not cancel. For floats
 * the error per edge is about 1e-16 |A-O| |B-O|: polygons far from
 * each other lose digits (use integer T).
 * Usage:
 *  typedef Point<ll> P;
 *  vector<vector<P>> polys = {{P(0,0),P(2,0),P(2,2),P(0,2)},
 *                             {P(1,1),P(3,1),P(3,3),P(1,3)}};
 *  double area = polyUnion(polys); // 7
 * Time: O(N^2), where N is the total number of points
 * Status: stress-tested, Submitted on ECNA 2017 Problem A
 */
#pragma once

#include "Point.h"
#include "sideOf.h"

template<class P> using PolyUnionT = GeoR<decltype(P().x)>;
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
	auto crs = [&](P a, P b) { // (a-O).cross(b-O)
		if constexpr (geoFloat<T>) { // products may cancel: geoDet
			typedef typename Geo<T>::W W;
			W x = W(a.x) - W(O.x), y = W(a.y) - W(O.y);
			return R(geoDet(x, W(b.y) - W(O.y), y, W(b.x) - W(O.x)));
		} else return R((a - O).cross(b - O));
	};
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
		ret += crs(A, B) * sum;
	}
	return ret / 2;
}
