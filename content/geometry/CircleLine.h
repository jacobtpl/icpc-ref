/**
 * Author: Victor Lecomte, chilli
 * Date: 2019-10-29
 * License: CC0
 * Source: https://vlecomte.github.io/cp-geo.pdf
 * Description: Finds the intersection between a circle and a line.
 * Returns a vector of either 0, 1, or 2 intersection points
 * (empty if a = b). P must be a floating point Point.
 * No epsilon is used. The count is exact for integer-valued input with
 * $|x|,|y|,r \le 10^9$ (x86 long double; the sign of
 * $r^2|ab|^2 - s^2$ is then computed exactly with fma).
 * Status: unit tested
 */

#pragma once

#include "Point.h"

template<class P>
vector<P> circleLine(P c, double r, P a, P b) {
	typedef long double ld; typedef decltype(a.x) T;
	static_assert(is_floating_point<T>::value, "use doubles");
	ld x = (ld)b.x - a.x, y = (ld)b.y - a.y, u = (ld)c.x - a.x,
	   v = (ld)c.y - a.y, L = x*x + y*y, s = x*v - y*u,
	   t = x*u + y*v, rr = (ld)r*r, w = s*s, N = rr*L - w;
	if (L == 0) return {};
	if (fabsl(N) < w * 1e-6L) // cancellation: exact sign (Kahan)
		N = fmal(rr, L, -w) + fmal(-s, s, w);
	if (N < 0) return {};
	auto at = [&](ld k) { return P(T(a.x+x*k), T(a.y+y*k)); };
	if (N == 0) return {at(t / L)};
	ld q = sqrtl(N);
	return {at((t - q) / L), at((t + q) / L)}; }
