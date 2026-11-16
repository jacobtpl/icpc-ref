/**
 * Author: chilli, Takanori MAEHARA
 * Date: 2019-10-31
 * License: CC0
 * Source: https://github.com/spaghetti-source/algorithm/blob/master/geometry/_geom.cc#L744
 * Description: Returns the area of the intersection of a circle with a
 * ccw polygon. Chords are built from the foot of the perpendicular from
 * the centre to each edge line, so small circles far from the vertices
 * are fine (e.g. r=1, coordinates 1e9).
 * P must be a floating point Point; r and the result have its
 * coordinate type. Everything (also the angles) is computed in
 * W = Geo<T>::W, which is T itself by default (see Point.h).
 * Usage:
 *  typedef Point<double> P;
 *  vector<P> sq = {P(0,0), P(2,0), P(2,2), P(0,2)};
 *  double a = circlePoly(P(0,0), 1.0, sq); // quarter disc: pi/4
 * Time: O(n)
 * Status: Tested on GNYR 2019 Gerrymandering, stress-tested
 */
#pragma once

#include "Point.h"

template<class P>
auto circlePoly(P c, decltype(P().x) r, const vector<P>& ps) {
	typedef decltype(c.x) T; typedef typename Geo<T>::W W;
	typedef Point<W> Q;
	static_assert(geoFloat<T>, "use a floating point Point");
	W R = r, r2 = R * R / 2;
	auto ang = [](Q p, Q q) {
		return gatan2(p.cross(q), p.dot(q)); };
	auto tri = [&](Q p, Q q) {
		Q d = q - p, u = p, v = q;
		W d2 = d.dist2(), cr = p.cross(d), h = gfabs(cr) / gsqrt(d2),
		  det = (R - h) * (R + h) / d2; // nan if p == q
		if (!(det > 0)) return ang(p, q) * r2;
		W a = d.dot(p) / d2, e = gsqrt(det),
		  s = -a - e, t = -a + e, k = -cr / d2;
		if (t < 0 || 1 <= s) return ang(p, q) * r2;
		// crossings = foot of the perpendicular d.perp()*k -+ d*e
		if (s > 0) u = Q(-d.y*k - d.x*e, d.x*k - d.y*e);
		if (t < 1) v = Q(-d.y*k + d.x*e, d.x*k + d.y*e);
		return ang(p,u) * r2 + u.cross(v)/2 + ang(v,q) * r2;
	};
	auto rel = [&](P p) { // p - c, formed in W
		return Q(W(p.x) - W(c.x), W(p.y) - W(c.y)); };
	W sum = 0;
	rep(i,0,sz(ps))
		sum += tri(rel(ps[i]), rel(ps[(i + 1) % sz(ps)]));
	return T(sum); }
