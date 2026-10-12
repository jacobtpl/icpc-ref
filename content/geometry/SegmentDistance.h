/**
 * Author: Ulf Lundstrom
 * Date: 2009-03-21
 * License: CC0
 * Source:
 * Description:\\
\begin{minipage}{65mm}
Returns the shortest distance between point p and the line segment from point s to e.
Works for any Point<T>; with Point<ll> (coordinates up to $10^9$) everything but the final division is exact, which double input cannot give near the segment.
\end{minipage}
\begin{minipage}{15mm}
\vspace{-10mm}
\includegraphics[width=\textwidth]{content/geometry/SegmentDistance}
\end{minipage}
 * Usage: 
 * 	Point<double> a, b(2,2), p(1,1);
 * 	bool onSegment = segDist(a,b,p) < 1e-10;
 * Status: tested
 */
#pragma once

#include "Point.h"

typedef Point<double> P; // default only; segDist is generic
template<class P> auto segDist(P s, P e, P p) {
	P v = e - s, w = p - s, u = p - e;
	typedef decltype(v.dist()) F;
	if (w.dot(v) <= 0) return w.dist(); // also s == e
	if (u.dot(v) >= 0) return u.dist();
	return abs((F)v.cross(w)) / v.dist();
}
