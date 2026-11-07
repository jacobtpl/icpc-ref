/**
 * Author: Jacob, Devin
 * Date: 2026-10-12
 * License: CC0
 * Description: Signed doubled area of $(a,c,X)$, where $X$
 * is the intersection of the nonparallel lines $ab$ and $cd$.
 * Integer coordinates must have absolute value at most $10^9$.
 * Computes $UV/D$, where $U=(c-a)\times(b-a)$,
 * $V=(c-a)\times(d-c)$ and $D=(b-a)\times(d-c)$.
 * Integer intermediates are exact; the result is approximate.
 * Use T=\texttt{\_\_float128} when later area sums may cancel.
 * Usage: auto area2 = lineInterArea2<__float128>(a,b,c,d);
 * Time: O(1)
 * Status: stress-tested against an independent rational oracle
 */
#pragma once

#include "Point.h"

template<class T = long double, class I>
T lineInterArea2(Point<I> a, Point<I> b,
		Point<I> c, Point<I> d) {
	static_assert(numeric_limits<I>::is_integer,
		"lineInterArea2 requires integer coordinates");
	static_assert(is_floating_point<T>::value ||
		is_same<T, __float128>::value,
		"lineInterArea2 requires a floating-point result");
	auto wide = [](Point<I> p) {
		return Point<__int128_t>(p.x, p.y);
	};
	auto A=wide(a), B=wide(b), C=wide(c), D=wide(d);
	auto u=A.cross(C,B), v=D.cross(A,C);
	auto den=(B-A).cross(D-C);
	assert(den != 0);
	return T(u*v) / T(den);
}
