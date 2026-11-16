/**
 * Author: Ulf Lundstrom
 * Date: 2009-04-08
 * License: CC0
 * Source:
 * Description: Returns the center of mass for a polygon of nonzero area
 *  (zero area or no points gives NaN, or throws for Fraction).
 *  Moments are taken relative to v[0], so the error is relative to the
 *  size of the polygon, not to its distance from the origin.
 *  T can be any type Point supports. The result is Point<GeoR<T>>
 *  (see Point.h): Point<double> for int and ll input, Point<T> for
 *  floating T and Fraction (exact). For another result type name it:
 *  polygonCenter<ll, long double>(v).
 *  Integer T: the sums are exact in W = Geo<T>::W and divided once at
 *  the end, valid while 4*D^3 fits in W for coordinate range D:
 *  ll: |x|,|y| < 2^30, int: < 5e5, __int128: < 1e12.
 * Usage:
 *  vector<Point<ll>> v = {Point<ll>(0,0), Point<ll>(4,0), Point<ll>(0,4)};
 *  Point<double> c = polygonCenter(v); // (4/3, 4/3)
 * Time: O(n)
 * Status: Tested
 */
#pragma once

#include "Point.h"

template<class T, class R = GeoR<T>>
Point<R> polygonCenter(const vector<Point<T>>& v) {
	typedef typename Geo<T>::W W;
	typedef conditional_t<geoInt<T>, R, W> C;
	W A = 0, X = 0, Y = 0;
	if (v.empty()) return Point<R>(R(C(A) / C(A)), R(C(A) / C(A)));
	rep(i,1,sz(v)-1) {
		W ax = W(v[i].x) - W(v[0].x), ay = W(v[i].y) - W(v[0].y),
		  bx = W(v[i+1].x) - W(v[0].x), by = W(v[i+1].y) - W(v[0].y),
		  c = ax * by - ay * bx;
		A += c; X += c * (ax + bx); Y += c * (ay + by);
	}
	C d = C(A * 3);
	return Point<R>(R(C(v[0].x) + C(X) / d), R(C(v[0].y) + C(Y) / d));
}
