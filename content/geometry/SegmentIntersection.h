/**
 * Author: Victor Lecomte, chilli
 * Date: 2019-04-27
 * License: CC0
 * Source: https://vlecomte.github.io/cp-geo.pdf
 * Description:\\
\begin{minipage}{65mm}
If a unique intersection point between the line segments going from s1 to e1 and from s2 to e2 exists then it is returned.
If no intersection point exists an empty vector is returned. If infinitely many exist a vector with 2 elements is returned, containing the endpoints of the common line segment.
For integer coordinates with absolute value at most $10^9$, classification is exact using 128-bit intermediates. The output defaults to Point<long double>; choose segInter<T> for another floating type or Fraction. Fraction output throws if a reduced coordinate does not fit. Rounding never changes the vector size, even if overlap endpoints round to equal points.
Floating-point input can misclassify nearly collinear segments. Fraction input is exact only while every intermediate fits. These inputs retain their scalar type. Requires C++17.
\end{minipage}
\begin{minipage}{15mm}
\includegraphics[width=\textwidth]{content/geometry/SegmentIntersection}
\end{minipage}
 * Usage:
 * auto inter = segInter(s1,e1,s2,e2);
 * if (sz(inter)==1)
 *   cout << "segments intersect at " << inter[0] << endl;
 * Status: stress-tested, tested on kattis:intersection
 */
#pragma once

#include "Point.h"
#include "OnSegment.h"
#include "../number-theory/Fraction.h"

template<class P> vector<P> segInter(P a, P b, P c, P d) {
	static_assert(!numeric_limits<decltype(a.x)>::is_integer,
		"Choose a floating scalar or Fraction for integer input");
	auto oa = c.cross(d, a), ob = c.cross(d, b),
	     oc = a.cross(b, c), od = a.cross(b, d);
	// Checks if intersection is single non-endpoint point.
	if (sgn(oa) * sgn(ob) < 0 && sgn(oc) * sgn(od) < 0)
		return {(a * ob - b * oa) / (ob - oa)};
	set<P> s;
	if (onSegment(c, d, a)) s.insert(a);
	if (onSegment(c, d, b)) s.insert(b);
	if (onSegment(a, b, c)) s.insert(c);
	if (onSegment(a, b, d)) s.insert(d);
	return {all(s)};
}

template<class T = long double, class I, enable_if_t<
	numeric_limits<I>::is_integer && (is_floating_point<T>::value ||
	is_same<T, __float128>::value || is_same<T, Fraction>::value), int> = 0>
vector<Point<T>> segInter(Point<I> a, Point<I> b, Point<I> c, Point<I> d) {
	using W = __int128_t;
	using P = Point<W>;
	auto wide = [](Point<I> p) { return P(p.x,p.y); };
	P A=wide(a), B=wide(b), C=wide(c), D=wide(d);
	auto ratio = [](W n, W den) -> T {
		if constexpr (is_same<T, Fraction>::value)
			return Fraction::fromWide(n,den);
		else return T(n)/T(den);
	};
	auto point = [&](P p, W den) {
		return Point<T>(ratio(p.x,den),ratio(p.y,den));
	};
	W oa=C.cross(D,A), ob=C.cross(D,B), oc=A.cross(B,C), od=A.cross(B,D);
	if (sgn(oa)*sgn(ob)<0 && sgn(oc)*sgn(od)<0)
		return {point(A*ob-B*oa,ob-oa)};
	set<P> s;
	if (onSegment(C,D,A)) s.insert(A);
	if (onSegment(C,D,B)) s.insert(B);
	if (onSegment(A,B,C)) s.insert(C);
	if (onSegment(A,B,D)) s.insert(D);
	vector<Point<T>> out;
	for (P p : s) out.push_back(point(p,1));
	return out;
}
