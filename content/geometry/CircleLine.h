/**
 * Author: Victor Lecomte, chilli
 * Date: 2019-10-29
 * License: CC0
 * Source: https://vlecomte.github.io/cp-geo.pdf
 * Description: Finds the intersection between a circle and a line.
 * Returns a vector of either 0, 1, or 2 intersection points
 * (empty if a = b). P must be a floating point Point; r has its
 * coordinate type. Everything is computed in W = Geo<T>::W (T itself by
 * default, see Point.h), so Point<long double> or Point<__float128>
 * give that precision throughout.
 * The sign of r^2|ab|^2 - s^2, which decides the count, is taken
 * with sgn (exact unless geoEps<W> is set) and is evaluated without
 * cancellation error. It is exact for integer-valued input while
 * r^2, |ab|^2 and s fit the mantissa of W: |x|,|y|,r <= 3e7 in
 * double, <= 1e9 in long double (as W or as T).
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
	  t = x*u + y*v, R = r, N = geoDet(R*R, L, s, s);
	if (L == 0 || sgn(N) < 0) return {};
	auto at = [&](W k) { return P(T(a.x+x*k), T(a.y+y*k)); };
	if (sgn(N) == 0) return {at(t / L)};
	W q = gsqrt(N);
	return {at((t - q) / L), at((t + q) / L)}; }
