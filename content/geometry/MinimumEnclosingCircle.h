/**
 * Author: Andrew He, chilli
 * Date: 2019-05-07
 * License: CC0
 * Source: folklore
 * Description: Computes the minimum circle that encloses a set of points.
 * Works in coordinates relative to one of the points, so the tolerance is
 * relative to the spread of the points, not to their distance from the
 * origin. Duplicate and collinear points are fine; empty input gives radius 0.
 * Points up to $10^{-12} r$ outside the circle are accepted.
 * Time: expected O(n)
 * Status: stress-tested
 */
#pragma once

#include "circumcircle.h"

pair<P, double> mec(vector<P> ps) {
	if (ps.empty()) return {P(), 0};
	static mt19937 rng((unsigned)time(0));
	shuffle(all(ps), rng);
	P s = ps[0], o;
	for (P& p : ps) p = p - s;
	double r = 0, EPS = 1 + 1e-12;
	rep(i,0,sz(ps)) if ((o - ps[i]).dist() > r * EPS) {
		o = ps[i], r = 0;
		rep(j,0,i) if ((o - ps[j]).dist() > r * EPS) {
			o = (ps[i] + ps[j]) / 2;
			r = (o - ps[i]).dist();
			rep(k,0,j) if ((o - ps[k]).dist() > r * EPS) {
				o = ccCenter(ps[i], ps[j], ps[k]);
				if (!isfinite(o.x + o.y)) { // collinear: farthest pair
					P a = ps[i] - ps[k], b = ps[j] - ps[k];
					o = ps[k] + (a.dist2() > b.dist2() ? a : b) / 2;
				}
				r = (o - ps[k]).dist();
			}
		}
	}
	return {o + s, r};
}
