/**
 * Author: Ulf Lundstrom
 * Date: 2009-04-11
 * License: CC0
 * Source: http://en.wikipedia.org/wiki/Circumcircle
 * Description: ccRadius = radius and ccCenter = center of the circle
 *  through points A, B, C (picture: circumcircle.pdf).
 *  Collinear points give inf/nan. Works on any floating Point and
 *  stays accurate for thin triangles; ccCenter also accepts
 *  Point<Fraction>. No eps is used.
 * Usage:
 *  typedef Point<double> P;
 *  P o = ccCenter(P(0,0), P(2,0), P(0,2));     // (1,1)
 *  double r = ccRadius(P(0,0), P(2,0), P(0,2)); // sqrt(2)
 * Status: tested
 */
#pragma once

#include "Point.h"

template<class P>
auto ccRadius(const P& A, const P& B, const P& C) {
	typedef decltype(A.x) T; typedef typename Geo<T>::W W;
	static_assert(geoFloat<T>, "use a floating point Point");
	W bx = W(C.x) - W(A.x), by = W(C.y) - W(A.y), cx = W(B.x) - W(A.x),
	  cy = W(B.y) - W(A.y), ex = bx - cx, ey = by - cy;
	return T(gsqrt((bx*bx + by*by) * (cx*cx + cy*cy) *
		(ex*ex + ey*ey)) / gfabs(geoDet(bx, cy, by, cx)) / 2); }
template<class P>
P ccCenter(const P& A, const P& B, const P& C) {
	typedef decltype(A.x) T; typedef typename Geo<T>::W W;
	if constexpr (geoFloat<T>) {
		W bx = W(C.x) - W(A.x), by = W(C.y) - W(A.y),
		  cx = W(B.x) - W(A.x), cy = W(B.y) - W(A.y),
		  bb = bx*bx + by*by, cc = cx*cx + cy*cy,
		  D = geoDet(bx, cy, by, cx) * 2;
		return P(T(A.x - geoDet(by, cc, cy, bb) / D),
		         T(A.y + geoDet(bx, cc, cx, bb) / D));
	} else {
		P b = C-A, c = B-A;
		return A + (b*c.dist2()-c*b.dist2()).perp()/b.cross(c)/2;
	} }
