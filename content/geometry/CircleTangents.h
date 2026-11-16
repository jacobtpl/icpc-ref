/**
 * Author: Victor Lecomte, chilli
 * Date: 2019-10-31
 * License: CC0
 * Source: https://vlecomte.github.io/cp-geo.pdf
 * Description: Finds the external tangents of two circles, or internal if r2 is negated.
 * Can return 0, 1, or 2 tangents -- 0 if one circle contains the other (or overlaps it, in the internal case, or if the circles are the same);
 * 1 if the circles are tangent to each other (in which case .first = .second and the tangent line is perpendicular to the line between the centers).
 * .first and .second give the tangency points at circle 1 and 2 respectively.
 * To find the tangents of a circle with a point set r2 to 0.
 * P is a floating point Point; the radii have its coordinate type.
 * The count is exact for integer-valued input up to 3e7 (double),
 * 1e9 (long double).
 * Eps: geoEps is compared to d^2 - (r1-r2)^2 (tangent circles).
 * Usage:
 *  typedef Point<double> P;
 *  auto t = tangents(P(0,0), 1.0, P(4,0), 1.0);   // 2 external
 *  auto in = tangents(P(0,0), 1.0, P(4,0), -1.0); // 2 internal
 * Status: tested
 */
#pragma once

#include "Point.h"

template<class P>
vector<pair<P, P>> tangents(P c1, decltype(P().x) r1,
		P c2, decltype(P().x) r2) {
	typedef decltype(c1.x) T; typedef typename Geo<T>::W W;
	static_assert(geoFloat<T>, "use a floating point Point");
	W x = W(c2.x) - W(c1.x), y = W(c2.y) - W(c1.y), R1 = r1, R2 = r2,
	  dr = R1 - R2, d2 = x*x + y*y, h2 = d2 - dr*dr;
	if (d2 == 0 || sgn(h2) < 0) return {};
	W h = h2 > 0 ? gsqrt(h2) : W(0);
	vector<pair<P, P>> out;
	for (int sign : {-1, 1}) {
		W vx = (x*dr - y*h*W(sign))/d2, vy = (y*dr + x*h*W(sign))/d2;
		out.push_back({P(T(c1.x + vx*R1), T(c1.y + vy*R1)),
		               P(T(c2.x + vx*R2), T(c2.y + vy*R2))});
	}
	if (sgn(h2) == 0) out.pop_back();
	return out; }
