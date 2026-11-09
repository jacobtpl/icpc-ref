/**
 * Author: Ulf Lundstrom
 * Date: 2009-03-21
 * License: CC0
 * Source: Basic math
 * Description:\\
\begin{minipage}{65mm}
Returns the signed distance between point p and the line containing points a and b. Positive value on left side and negative on right as seen from a towards b. a==b gives nan. P is supposed to be Point<T> or Point3D<T> where T is e.g. double or long long. The cross product is evaluated in T: exact for long long coordinates up to $10^9$ (prefer that to double for integer input, a double cross product cancels for points near the line). The result has the floating type of P (double for integers). Using Point3D will always give a non-negative distance. For Point3D, call .dist on the result of the cross product.
\end{minipage}
\begin{minipage}{15mm}
\includegraphics[width=\textwidth]{content/geometry/lineDistance}
\end{minipage}
 * Status: tested
 */
#pragma once

#include "Point.h"

template<class P>
auto lineDist(const P& a, const P& b, const P& p) {
	typedef decltype((b-a).dist()) F;
	return (F)(b-a).cross(p-a)/(b-a).dist();
}
