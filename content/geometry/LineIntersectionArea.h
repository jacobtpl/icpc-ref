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
 * The result type defaults to Geo<I>::F (see Point.h; double for
 * int and long long); choose another with lineInterArea2<T>.
 * Use T = long double or __float128 when later area sums may cancel.
 * Usage: auto area2 = lineInterArea2<__float128>(a,b,c,d);
 * Time: O(1)
 * Status: stress-tested against an independent rational oracle
 */
#pragma once

#include "Point.h"

template<class R = void, class I>
auto lineInterArea2(Point<I> a, Point<I> b,
		Point<I> c, Point<I> d) {
	typedef conditional_t<is_void_v<R>, typename Geo<I>::F, R> T;
	static_assert(geoInt<I>,
		"lineInterArea2 requires integer coordinates");
	static_assert(geoFloat<T>,
		"lineInterArea2 requires a floating-point result");
	auto wide = [](Point<I> p) { // U*V has degree 4: beyond Geo<I>::W
		return Point<__int128>(p.x, p.y);
	};
	auto A=wide(a), B=wide(b), C=wide(c), D=wide(d);
	auto u=A.cross(C,B), v=D.cross(A,C);
	auto den=(B-A).cross(D-C);
	assert(den != 0);
	return T(u*v) / T(den);
}
