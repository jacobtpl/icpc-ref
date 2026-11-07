/**
 * Author: Ulf Lundstrom
 * Date: 2009-03-21
 * License: CC0
 * Source: tinyKACTL
 * Description: Returns twice the signed area of a polygon.
 *  Clockwise enumeration gives negative area. Empty input gives 0.
 *  Coordinates are taken relative to v[0], so a small polygon far from
 *  the origin does not cancel in floating point.
 *  With T = ll every term is exact for $|x|,|y| < 2^{30}$ (about $1.07\cdot 10^9$)
 *  and so is the result; the running sum of a non-convex polygon can
 *  exceed the result, so near that limit accumulate in a wider type:
 *  \texttt{polygonArea2<ll, \_\_int128>(v)}.
 * Status: Stress-tested and tested on kattis:polygonarea
 */
#pragma once

#include "Point.h"

template<class T, class A = T>
A polygonArea2(const vector<Point<T>>& v) {
	A a = 0;
	rep(i,1,sz(v)-1) a += (A)(v[i]-v[0]).cross(v[i+1]-v[0]);
	return a;
}
