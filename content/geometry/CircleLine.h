/**
 * Author: Victor Lecomte, chilli
 * Date: 2019-10-29
 * License: CC0
 * Source: https://vlecomte.github.io/cp-geo.pdf
 * Description: Finds the intersection between a circle and a line.
 * Returns a vector of either 0, 1, or 2 intersection points
 * (empty if a = b). P is a floating point Point; r has its
 * coordinate type. The count is exact for integer-valued input up to
 * 3e7 (double), 1e9 (long double). Points are built from the foot
 * of the perpendicular through c, so their error is relative to r,
 * also when a and b are far away.
 * Eps: geoEps is compared to r^2|ab|^2 - s^2 with s = (b-a)x(c-a),
 * a length^4, so scale it by |ab|^2 (or normalise b-a).
 * Usage:
 *  typedef Point<double> P;
 *  vector<P> v = circleLine(P(0,0), 5.0, P(-10,3), P(10,3));
 *  // v = {(-4,3), (4,3)}, ordered along a -> b
 * Status: unit tested
 */

#pragma once

#include "Point.h"

template<class P>
vector<P> circleLine(P c, decltype(P().x) r, P a, P b) {
	typedef decltype(a.x) T; typedef typename Geo<T>::W W;
	static_assert(geoFloat<T>, "use a floating point Point");
	W x = W(b.x) - W(a.x), y = W(b.y) - W(a.y), u = W(c.x) - W(a.x),
	  v = W(c.y) - W(a.y), L = x*x + y*y, s = x*v - y*u,
	  R = r, N = geoDet(R*R, L, s, s);
	if (L == 0 || sgn(N) < 0) return {};
	W fx = c.x + y*s/L, fy = c.y - x*s/L; // foot of c on the line
	auto at = [&](W k) { return P(T(fx+x*k), T(fy+y*k)); };
	if (sgn(N) == 0) return {at(0)};
	W q = gsqrt(N) / L;
	return {at(-q), at(q)}; }
