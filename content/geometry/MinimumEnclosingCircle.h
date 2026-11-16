/**
 * Author: Andrew He, chilli
 * Date: 2019-05-07
 * License: CC0
 * Source: folklore
 * Description: Computes the minimum circle that encloses a set of points.
 * P must be a floating point Point; returns {center, radius}.
 * Works in coordinates relative to one of the points, so the tolerance is
 * relative to the spread of the points, not to their distance from the
 * origin. Duplicate and collinear points are fine; empty input gives radius 0.
 * Points up to eps * r outside the circle are accepted. eps is
 * relative (it does not depend on the size of the coordinates) and
 * only has to cover rounding in the circle through 3 points:
 * 1e-12 suits double, 1e-15 long double, 1e-28 __float128. Too small
 * costs time (extra rebuilds), not correctness; too large returns a
 * circle up to that much too small.
 * Usage:
 *  typedef Point<double> P;
 *  vector<P> v = {P(0,0), P(2,0), P(1,1)};
 *  auto [o, r] = mec(v);                  // o = (1,0), r = 1
 *  auto res = mec(vector<Point<long double>>{...}, 1e-15L);
 * Time: expected O(n)
 * Status: stress-tested
 */
#pragma once

#include "circumcircle.h"

template<class P>
pair<P, decltype(P().x)> mec(vector<P> ps,
		decltype(P().x) eps = 1e-12) {
	typedef decltype(P().x) T;
	static_assert(geoFloat<T>, "use a floating point Point");
	if (ps.empty()) return {P(), 0};
	static mt19937 rng((unsigned)time(0));
	shuffle(all(ps), rng);
	P s = ps[0], o;
	for (P& p : ps) p = p - s;
	T r = 0, EPS = 1 + eps;
	rep(i,0,sz(ps)) if ((o - ps[i]).dist() > r * EPS) {
		o = ps[i], r = 0;
		rep(j,0,i) if ((o - ps[j]).dist() > r * EPS) {
			o = (ps[i] + ps[j]) / 2;
			r = (o - ps[i]).dist();
			rep(k,0,j) if ((o - ps[k]).dist() > r * EPS) {
				o = ccCenter(ps[i], ps[j], ps[k]);
				T f = o.x + o.y;
				if (!(f - f == 0)) { // inf/nan: collinear, farthest pair
					P a = ps[i] - ps[k], b = ps[j] - ps[k];
					o = ps[k] + (a.dist2() > b.dist2() ? a : b) / 2;
				}
				r = (o - ps[k]).dist();
			}
		}
	}
	return {o + s, r};
}
