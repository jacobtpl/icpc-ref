/**
 * Author: Ulf Lundstrom
 * Date: 2009-04-08
 * License: CC0
 * Source:
 * Description: Returns the center of mass for a polygon of nonzero area
 *  (zero area gives NaN, or throws for Fraction). Moments are taken
 *  relative to v[0], so the error is relative to the size of the polygon,
 *  not to its distance from the origin; Fraction input is exact. The ll
 *  overload sums exactly in \_\_int128 for $|x|,|y| < 2^{30}$ and divides once.
 * Time: O(n)
 * Status: Tested
 */
#pragma once

#include "Point.h"

typedef Point<double> P;
template<class P>
P polygonCenter(const vector<P>& v) {
	typedef decltype(P().x) T;
	T A = 0;
	if (v.empty()) return P(A / A, A / A);
	P res(0, 0);
	rep(i,1,sz(v)-1) {
		P a = v[i] - v[0], b = v[i+1] - v[0];
		T c = a.cross(b);
		res = res + (a + b) * c; A += c;
	}
	return v[0] + res / A / 3;
}
Point<long double> polygonCenter(const vector<Point<ll>>& v) {
	typedef long double ld;
	if (v.empty()) return Point<ld>(NAN, NAN);
	__int128 A = 0, X = 0, Y = 0;
	rep(i,1,sz(v)-1) {
		Point<ll> a = v[i] - v[0], b = v[i+1] - v[0];
		__int128 c = a.cross(b);
		A += c; X += c * (a.x + b.x); Y += c * (a.y + b.y);
	}
	return Point<ld>((ld)v[0].x + (ld)X / (ld)(3 * A),
		(ld)v[0].y + (ld)Y / (ld)(3 * A));
}
