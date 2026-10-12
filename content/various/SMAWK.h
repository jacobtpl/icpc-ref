/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: Aggarwal, Klawe, Moran, Shor, Wilber (1987), "Geometric applications of a matrix-searching algorithm"; own implementation
 * Description: Leftmost row minima of an implicit totally monotone $n \times m$ matrix $A[i][j] = f(i, j)$:
 *  for all $i < i'$, $j < j'$, $A[i][j] > A[i][j']$ must imply $A[i'][j] > A[i'][j']$
 *  (so the argmins are non-decreasing). This holds for Monge matrices:
 *  $A[i][j] + A[i'][j'] \le A[i][j'] + A[i'][j]$.
 *  Returns the 0-indexed column of the leftmost minimum of each row. Requires $m \ge 1$ if $n \ge 1$.
 *  Values are only compared with $<$ and $>$. Negating $f$ gives maxima only if $-A$
 *  satisfies the condition (e.g. inverse Monge $A$); for row maxima of a Monge $A$ use
 *  $g(i, j) = -f(i, m-1-j)$: row $i$'s rightmost maximum is at column $m-1-$ans$[i]$.
 *  Fewer calls to $f$ than D\&C, but rarely faster unless $f$ is expensive.
 * Time: O(N + M) calls to $f$
 * Usage: vi opt = smawk(n, m, [\&](int i, int j) { return cost(i, j); });
 * Status: stress-tested
 */
#pragma once

template<class F>
void smawkRec(F& f, vi& ans, int s, vi c) {
	int n = sz(ans) / s, j = 0; // rows s-1, 2s-1, ...
	if (!n) return;
	vi d;
	for (int x : c) {
		while (!d.empty() &&
			f(sz(d)*s-1, d.back()) > f(sz(d)*s-1, x))
			d.pop_back();
		if (sz(d) < n) d.push_back(x);
	}
	smawkRec(f, ans, 2 * s, d);
	for (int k = 0; k < n; k += 2) {
		int r = (k+1)*s-1, hi = k+1 < n ? ans[r+s] : d.back();
		ans[r] = d[j];
		while (d[j] != hi)
			if (f(r, d[++j]) < f(r, ans[r])) ans[r] = d[j];
	}
}
template<class F> vi smawk(int n, int m, F f) {
	vi ans(n), c(m);
	iota(all(c), 0);
	smawkRec(f, ans, 1, c);
	return ans;
}
