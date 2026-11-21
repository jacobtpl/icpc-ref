/**
 * Author: Ulf Lundstrom
 * Date: 2009-03-21
 * License: CC0
 * Source:
 * Description:
Returns the shortest distance between point p and the line segment from point s to e.
Works for any Point<T>; with Point<ll> (coordinates up to 1e9) everything but the final division is exact, which double input cannot give near the segment.
 * Usage: 
 * 	Point<double> a, b(2,2), p(1,1);
 * 	bool onSegment = segDist(a,b,p) < 1e-10;
 * Status: tested
 */
#pragma once

#include "Point.h"

template<class P> auto segDist(P s, P e, P p) {
	P v = e - s, w = p - s, u = p - e;
	typedef decltype(v.dist()) F;
	if (w.dot(v) <= 0) return w.dist(); // also s == e
	if (u.dot(v) >= 0) return u.dist();
	return gfabs((F)v.cross(w)) / v.dist();
}
