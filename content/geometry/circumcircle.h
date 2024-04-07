/**
 * Author: Ulf Lundstrom
 * Date: 2009-04-11
 * License: CC0
 * Source: http://en.wikipedia.org/wiki/Circumcircle
 * Description:\\
\begin{minipage}{75mm}
ccRadius = radius and ccCenter = center of circle through points A, B, C.
\end{minipage}
\begin{minipage}{15mm}
\vspace{-2mm}
\includegraphics[width=\textwidth]{content/geometry/circumcircle}
\end{minipage}
 * Status: tested
 */
#pragma once

#include "Point.h"

typedef Point<double> P;
double ccRadius(const P& A, const P& B, const P& C) {
	return (B-A).dist()*(C-B).dist()*(A-C).dist()/
			abs((B-A).cross(C-A))/2; }
P ccCenter(const P& A, const P& B, const P& C) {
	P b = C-A, c = B-A;
	return A + (b*c.dist2()-c*b.dist2()).perp()/b.cross(c)/2; }
