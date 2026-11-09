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
 * P must be a floating point Point. No epsilon is used; the count is exact for
 * integer-valued input with $|x|,|y|,r \le 10^9$ (x86 long double).
 * Status: tested
 */
#pragma once

#include "Point.h"

template<class P>
vector<pair<P, P>> tangents(P c1, double r1, P c2, double r2) {
	typedef long double ld; typedef decltype(c1.x) T;
	static_assert(is_floating_point<T>::value, "use doubles");
	ld x = (ld)c2.x - c1.x, y = (ld)c2.y - c1.y, dr = (ld)r1 - r2,
	   d2 = x*x + y*y, h2 = d2 - dr*dr;
	if (d2 == 0 || h2 < 0)  return {};
	ld h = sqrtl(h2);
	vector<pair<P, P>> out;
	for (ld sign : {-1, 1}) {
		ld vx = (x*dr - y*h*sign)/d2, vy = (y*dr + x*h*sign)/d2;
		out.push_back({P(T(c1.x + vx*r1), T(c1.y + vy*r1)),
		               P(T(c2.x + vx*r2), T(c2.y + vy*r2))});
	}
	if (h2 == 0) out.pop_back();
	return out; }
