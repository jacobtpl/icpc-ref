/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: own work; standard sweep line over a set of segments ordered by y
 * Description: Offline point location. For each query point returns the index
 * of the first segment hit by a ray going straight down from it (starting at
 * the point itself), or $-1$. A query lying on segments returns any of them.
 * Segments must have nonzero length and may only meet at a point that is an
 * endpoint of one of them (identical duplicates are also fine).
 * Ties in $x$ are broken by $y$, i.e. the ray leans infinitesimally towards
 * $+x$: with endpoints $a < b$ (as Points), a segment can only be hit if
 * $a \le q \le b$. So a segment whose right end is straight below $q$ is not
 * hit, one whose left end is there is hit, and a vertical segment is hit only
 * by queries on it. A query on no segment lies in the face just above the
 * returned segment. Coordinates must be at most $10^9$ in absolute value.
 * Time: O((n + q) \log (n + q))
 * Usage: vi below = pointLocation(segs, queries);
 * Status: stress-tested
 */
#pragma once

#include "Point.h"

typedef Point<ll> P;
vi pointLocation(vector<pair<P, P>> s, const vector<P>& q) {
	int n = sz(s), m = sz(q), e = -1;
	for (auto& [a, b] : s) if (b < a) swap(a, b);
	for (P p : q) s.emplace_back(p, p);
	auto cmp = [&](int i, int j) { // is i below j?
		auto [a, b] = s[i]; auto [c, d] = s[j];
		if (a < c) {
			ll x = a.cross(b, c);
			return x ? x > 0 : a.cross(b, d) > 0;
		}
		ll x = c.cross(d, a);
		return x ? x < 0 : c.cross(d, b) < 0;
	};
	set<int, decltype(cmp)> act(cmp);
	vector<tuple<P, int, int>> ev; // 0 = end, 1 = start, 2 = query
	rep(i,0,n) ev.emplace_back(s[i].second, 0, i),
		ev.emplace_back(s[i].first, 1, i);
	rep(i,n,n+m) ev.emplace_back(s[i].first, 2, i);
	sort(all(ev));
	vi res(m);
	for (auto [p, t, i] : ev) {
		if (t == 0) act.erase(i), e = i;
		else if (t == 1) act.insert(i);
		else {
			auto it = act.upper_bound(i);
			res[i-n] = e >= 0 && s[e].second == p ? e :
				it == act.begin() ? -1 : *--it;
		}
	}
	return res;
}
