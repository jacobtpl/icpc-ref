/**
 * Author: chilli, Takanori MAEHARA
 * Date: 2019-10-31
 * License: CC0
 * Source: https://github.com/spaghetti-source/algorithm/blob/master/geometry/_geom.cc#L744
 * Description: Returns the area of the intersection of a circle with a
 * ccw polygon. Chords are built from the foot of the perpendicular from
 * the centre to each edge line, so small circles far from the vertices
 * are fine (e.g. $r=1$, coordinates $10^9$).
 * Time: O(n)
 * Status: Tested on GNYR 2019 Gerrymandering, stress-tested
 */
#pragma once

#include "../../content/geometry/Point.h"

typedef Point<double> P;
double circlePoly(P c, double r, vector<P> ps) {
	typedef long double ld;
	auto ang = [](P p, P q) {
		return atan2(p.cross(q), p.dot(q)); };
	auto tri = [&](P p, P q) {
		auto r2 = r * r / 2;
		P d = q - p, u = p, v = q;
		ld d2 = (ld)d.x*d.x + (ld)d.y*d.y,
		   cr = (ld)p.x*d.y - (ld)p.y*d.x, h = fabsl(cr) / sqrtl(d2),
		   det = (r - h) * (r + h) / d2; // nan if p == q
		if (!(det > 0)) return ang(p, q) * r2;
		ld a = ((ld)d.x*p.x + (ld)d.y*p.y) / d2, e = sqrtl(det),
		   s = -a - e, t = -a + e, k = -cr / d2;
		if (t < 0 || 1 <= s) return ang(p, q) * r2;
		// crossings = foot of the perpendicular d.perp()*k -+ d*e
		if (s > 0) u = P(double(-d.y*k-d.x*e), double(d.x*k-d.y*e));
		if (t < 1) v = P(double(-d.y*k+d.x*e), double(d.x*k+d.y*e));
		return ang(p,u) * r2 + u.cross(v)/2 + ang(v,q) * r2;
	};
	auto sum = 0.0;
	rep(i,0,sz(ps))
		sum += tri(ps[i] - c, ps[(i + 1) % sz(ps)] - c);
	return sum; }
