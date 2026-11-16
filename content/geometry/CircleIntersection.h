/**
 * Author: Simon Lindholm
 * Date: 2015-09-01
 * License: CC0
 * Description: Computes the pair of points at which two circles intersect. Returns false in case of no intersection.
 * Asserts on identical circles. Tangent circles give two equal points.
 * P must be a floating point Point; the radii have its coordinate type.
 * Everything is computed in W = Geo<T>::W (T itself by default, see
 * Point.h). The yes/no answer is decided with sgn (exact unless
 * geoEps<W> is set); it is exact for integer-valued input while squared
 * distances fit the mantissa of W: |x|,|y|,r <= 3e7 in double,
 * <= 1e9 in long double.
 * Usage:
 *  typedef Point<double> P;
 *  pair<P, P> out;
 *  if (circleInter(P(0,0), P(4,0), 3.0, 3.0, &out))
 *    cout << out.first << " " << out.second << endl; // (2,-2.236) (2,2.236)
 * Status: stress-tested
 */
#pragma once

#include "Point.h"

template<class P>
bool circleInter(P a, P b, decltype(P().x) r1,
		decltype(P().x) r2, pair<P, P>* out) {
	typedef decltype(a.x) T; typedef typename Geo<T>::W W;
	static_assert(geoFloat<T>, "use a floating point Point");
	if (a == b) { assert(r1 != r2); return false; }
	W x = W(b.x) - W(a.x), y = W(b.y) - W(a.y), d2 = x*x + y*y,
	  sum = W(r1) + W(r2), dif = W(r1) - W(r2),
	  u = sum*sum - d2, v = d2 - dif*dif; // h2 = u*v/(4*d2)
	if (sgn(u) < 0 || sgn(v) < 0) return false;
	W p = (1 + sum*dif/d2) / 2, uv = u*v,
	  h = (uv > 0 ? gsqrt(uv) : W(0)) / (d2*2),
	  mx = a.x + x*p, my = a.y + y*p;
	*out = {P(T(mx - y*h), T(my + x*h)),
	        P(T(mx + y*h), T(my - x*h))};
	return true; }
