/**
 * Author: Ulf Lundstrom
 * Date: 2009-04-11
 * License: CC0
 * Source: http://en.wikipedia.org/wiki/Circumcircle
 * Description:\\
\begin{minipage}{65mm}
ccRadius = radius and ccCenter = center of circle through points A, B, C.
\end{minipage}
\begin{minipage}{15mm}
\vspace{-2mm}
\includegraphics[width=\textwidth]{content/geometry/circumcircle}
\end{minipage}
 Collinear points give inf/nan. Differences and products are taken in long
 double with fma on cancellation, so thin triangles with integer coordinates
 up to $10^9$ are accurate.
 * Status: tested
 */
#pragma once

#include "Point.h"

typedef Point<double> P;
// a*b - c*d, accurate also when the products cancel (Kahan)
long double ccDet(long double a, long double b,
		long double c, long double d) {
	long double w = c*d, x = a*b - w;
	if (fabsl(x) > fabsl(w) * 1e-3L) return x;
	return fmal(a, b, -w) + fmal(-c, d, w); }
double ccRadius(const P& A, const P& B, const P& C) {
	typedef long double ld;
	ld bx = (ld)C.x - A.x, by = (ld)C.y - A.y, cx = (ld)B.x - A.x,
	   cy = (ld)B.y - A.y, ex = bx - cx, ey = by - cy;
	return double(sqrtl((bx*bx + by*by) * (cx*cx + cy*cy) *
		(ex*ex + ey*ey)) / fabsl(ccDet(bx, cy, by, cx)) / 2); }
P ccCenter(const P& A, const P& B, const P& C) {
	typedef long double ld;
	ld bx = (ld)C.x - A.x, by = (ld)C.y - A.y, cx = (ld)B.x - A.x,
	   cy = (ld)B.y - A.y, bb = bx*bx + by*by, cc = cx*cx + cy*cy,
	   D = ccDet(bx, cy, by, cx) * 2;
	return P(double(A.x - ccDet(by, cc, cy, bb) / D),
	         double(A.y + ccDet(bx, cc, cx, bb) / D)); }
