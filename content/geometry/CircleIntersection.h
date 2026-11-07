/**
 * Author: Simon Lindholm
 * Date: 2015-09-01
 * License: CC0
 * Description: Computes the pair of points at which two circles intersect. Returns false in case of no intersection.
 * Asserts on identical circles. No epsilon is used: tangent circles give two equal points.
 * The yes/no answer is exact for integer-valued input with $|x|,|y|,r \le 10^9$
 * (squared distances then fit the 64-bit mantissa of x86 long double).
 * Status: stress-tested
 */
#pragma once

#include "Point.h"

typedef Point<double> P;
bool circleInter(P a,P b,double r1,double r2,pair<P, P>* out) {
	typedef long double ld;
	if (a == b) { assert(r1 != r2); return false; }
	ld x = (ld)b.x - a.x, y = (ld)b.y - a.y, d2 = x*x + y*y,
	   sum = (ld)r1 + r2, dif = (ld)r1 - r2,
	   u = sum*sum - d2, v = d2 - dif*dif; // h2 = u*v/(4*d2)
	if (u < 0 || v < 0) return false;
	ld p = (1 + sum*dif/d2) / 2, h = sqrtl(u*v) / (d2*2),
	   mx = a.x + x*p, my = a.y + y*p;
	*out = {P(double(mx - y*h), double(my + x*h)),
	        P(double(mx + y*h), double(my - x*h))};
	return true; }
