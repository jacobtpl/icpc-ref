/**
 * Author: Ulf Lundstrom
 * Date: 2009-03-21
 * License: CC0
 * Source:
 * Description: Returns a vector with the vertices of a polygon with
 *  everything to the left of the line going from s to e cut away
 *  (picture: PolygonCut.pdf).
 *  T can be any type Point supports. Side tests are evaluated in T
 *  (exact for integers and Fraction; watch for overflow) and compared
 *  exactly. The result is vector<Point<GeoR<T>>> (see Point.h):
 *  Point<double> for int and ll input, so new vertices are never
 *  truncated; Point<T> for floating T and Fraction.
 * Usage:
 * 	typedef Point<double> P;
 * 	vector<P> p = ...;
 * 	p = polygonCut(p, P(0,0), P(1,0)); // keeps the part with y <= 0
 * Status: tested but not extensively
 */
#pragma once

#include "Point.h"

template<class P>
auto polygonCut(const vector<P>& poly, P s, P e) {
	typedef decltype(s.x) T;
	typedef GeoC<T> C; typedef Point<GeoR<T>> Q;
	vector<Q> res;
	rep(i,0,sz(poly)) {
		Point<T> cur = poly[i], prev = i ? poly[i-1] : poly.back();
		T a = s.cross(e, cur), b = s.cross(e, prev);
		if ((a < 0) != (b < 0)) {
			C k = C(a) / (C(a) - C(b)), x = C(cur.x), y = C(cur.y);
			res.push_back(Q(GeoR<T>(x + (C(prev.x) - x) * k),
			                GeoR<T>(y + (C(prev.y) - y) * k)));
		}
		if (a < 0) res.push_back(Q(GeoR<T>(cur.x), GeoR<T>(cur.y)));
	}
	return res;
}
