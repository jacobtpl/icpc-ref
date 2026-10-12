/**
 * Author: Victor Lecomte, chilli
 * Date: 2019-05-05
 * License: CC0
 * Source: https://vlecomte.github.io/cp-geo.pdf
 * Description:\\
\begin{minipage}{65mm}
If a unique intersection point of the lines going through s1,e1 and s2,e2 exists \{1, point\} is returned.
If no intersection point exists \{0, (0,0)\} is returned and if infinitely many exists \{-1, (0,0)\} is returned.
A degenerate line (s1 = e1) is treated as a point: -1 if it lies on the other line, else 0 (never 1).
The point is built as an offset from s1, so doubles lose no digits to the absolute position.
For Point<ll> coordinates must be at most $10^9$ in absolute value and the answer must fit in ll (triple products are done in \_\_int128); each coordinate of a non-integer intersection is truncated toward s1.
With doubles the parallel test is exact, so nearly parallel lines give a far away point.
\end{minipage}
\begin{minipage}{15mm}
\includegraphics[width=\textwidth]{content/geometry/lineIntersection}
\end{minipage}
 * Usage:
 * 	auto res = lineInter(s1,e1,s2,e2);
 * 	if (res.first == 1)
 * 		cout << "intersection point at " << res.second << endl;
 * Status: stress-tested, and tested through half-plane tests
 */
#pragma once

#include "Point.h"

template<class P>
pair<int, P> lineInter(P s1, P e1, P s2, P e2) {
	auto d = (e1 - s1).cross(e2 - s2);
	if (d == 0) { // parallel, or a line is a single point
		if (s1 == e1) swap(s1, s2), swap(e1, e2);
		return {-(s1 == e1 ? s1 == s2 :
			s1.cross(e1, s2) == 0), P(0, 0)};
	}
	auto q = s2.cross(e2, s1); P v = e1 - s1;
	return {1, s1 + P(mulDiv(v.x, q, d), mulDiv(v.y, q, d))};
}
