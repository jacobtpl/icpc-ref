/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: own work; idea from IOI 2016 Aliens, https://codeforces.com/blog/entry/49691
 * Description: Lambda optimisation (Aliens trick). Let $c(j)$ be the min cost using
 *  exactly $j$ items, convex in $j$ (integer values). $f(\lambda)$ must return
 *  $(\min_j c(j) + \lambda j, \text{ any optimal } j)$, i.e. solve the problem with
 *  no limit on the number of items, each item costing $\lambda$ extra.
 *  Returns $c(k)$; $k$ must be feasible.
 *  $[lo, hi]$ must contain a $\lambda$ for which $k$ is optimal; $\pm\max_j |c(j)-c(j-1)|$
 *  always works. $\lambda k$ and the values must fit in ll.
 *  Ties are handled: $f$ may break them arbitrarily and need never return count $k$.
 *  To maximise a concave $c$, negate the costs and the result.
 * Usage:
	aliens(-C, C, k, [\&](ll lam) {
		... return pair<ll, int>(val, cnt); });
 * Time: O(\log(hi-lo)) calls to $f$
 * Status: stress-tested
 */
#pragma once

template<class F>
ll aliens(ll lo, ll hi, ll k, F f) {
	ll res = LLONG_MIN;
	while (lo <= hi) {
		ll m = lo + (hi - lo) / 2;
		auto [v, c] = f(m);
		res = max(res, v - m * k);
		if (c > k) lo = m + 1;
		else hi = m - 1;
	}
	return res;
}
