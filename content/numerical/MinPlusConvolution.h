/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: own work; the folklore merge of slopes (Minkowski sum) and
 * divide and conquer over monotone row minima (cf. DivideAndConquerDP.h)
 * Description: Min-plus convolution $c[k] = \min_{i+j=k} a[i]+b[j]$,
 * $0 \le k < |a|+|b|-1$ (empty if either input is empty).
 * $v$ is convex if $v[i+1]-v[i]$ is non-decreasing.
 * \texttt{minPlusCC} needs both $a$ and $b$ convex, \texttt{minPlus} needs
 * only $a$ convex ($b$ arbitrary). Every sum $a[i]+b[j]$ must fit in a ll;
 * nothing else is computed. For max-plus with concave input, negate everything.
 * Time: O(N) for \texttt{minPlusCC}, O(N \log N) for \texttt{minPlus}, where $N = |a|+|b|$
 * Status: stress-tested
 */
#pragma once

typedef vector<ll> vl;
vl minPlusCC(const vl& a, const vl& b) {
	int n = sz(a), m = sz(b), i = 0, j = 0;
	if (!n || !m) return {};
	vl c(n + m - 1);
	rep(k,0,n+m-1) {
		c[k] = a[i] + b[j];
		if (j == m-1 || (i < n-1 && a[i+1] + b[j] < a[i] + b[j+1]))
			i++;
		else j++;
	}
	return c;
}

vl minPlus(const vl& a, const vl& b) {
	int n = sz(a), m = sz(b);
	if (!n || !m) return {};
	vl c(n + m - 1);
	auto rec = [&](auto& f, int l, int r, int lo, int hi) -> void {
		if (l >= r) return;
		int k = (l + r) / 2, o = -1;
		rep(j, max(lo, k-n+1), min(hi, k) + 1) {
			ll v = a[k-j] + b[j];
			if (o < 0 || v < c[k]) c[k] = v, o = j;
		}
		f(f, l, k, lo, o), f(f, k+1, r, o, hi);
	};
	rec(rec, 0, n+m-1, 0, m-1);
	return c;
}
