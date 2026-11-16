/**
 * Author: Simon Lindholm
 * Date: 2019-04-17
 * License: CC0
 * Source: https://codeforces.com/blog/entry/58747
 * Description: Finds the closest pair of points (at least 2 points).
 * Works on any Point<T>. Squared distances are evaluated in T, so
 * 2 * (coordinate range)^2 must fit: |x|,|y| <= 1e9 for long long,
 * <= 1.6e4 for int. The strip width uses a square root in Geo<T>::F
 * (see Point.h), rounded up, so integer input stays exact.
 * Usage:
 *  typedef Point<ll> P;
 *  vector<P> v = {P(0,0), P(5,5), P(1,1), P(9,0)};
 *  pair<P, P> best = closest(v); // (0,0) and (1,1)
 * Time: O(n log n)
 * Status: stress-tested
 */
#pragma once

#include "Point.h"

template<class P>
pair<P, P> closest(vector<P> v) {
	typedef decltype(P().x) T; typedef typename P::F F;
	assert(sz(v) > 1);
	set<P> S;
	sort(all(v), [](P a, P b) { return a.y < b.y; });
	pair<T, pair<P, P>> ret{(v[0] - v[1]).dist2(), {v[0], v[1]}};
	int j = 0, i = 0;
	for (P p : v) {
		P d{T(gsqrt(F(ret.first))) + !geoFloat<T>, 0};
		// too far: tested on the differences dist2 uses, so rounding
		// of p.y - d or p.x - d can not drop a closer point
		auto far = [&](T a, T b) { return (a-b)*(a-b) > ret.first; };
		while (j < i && far(p.y, v[j].y)) S.erase(v[j++]);
		i++;
		auto lo = S.lower_bound(p - d), hi = S.upper_bound(p + d);
		while (lo != S.begin() && !far(p.x, prev(lo)->x)) --lo;
		while (hi != S.end() && !far(p.x, hi->x)) ++hi;
		for (; lo != hi; ++lo)
			ret = min(ret, {(*lo - p).dist2(), {*lo, p}});
		S.insert(p);
	}
	return ret.second;
}
