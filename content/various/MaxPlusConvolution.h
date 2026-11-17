/**
 * Author: Devin
 * Date: 2026-11-17
 * License: CC0
 * Source: folklore (Minkowski sum of concave sequences; monotone argmax);
 *  UCup 5 Stage 3 K
 * Description: c[k] = max over i+j=k of a[i] + b[j], sz(c) = sz(a)+sz(b)-1.
 *  Concave means a[i+1]-a[i] is non-increasing (e.g. best value using
 *  exactly i items, when taking items greedily is optimal).
 *  maxPlusCC: both concave, result is concave.
 *  maxPlusCA: a concave, b arbitrary.
 *  Both inputs must be non-empty. For min-plus with convex inputs, negate
 *  inputs and output.
 * Time: O(N+M) for maxPlusCC, O((N+M) log(N+M)) for maxPlusCA
 * Status: stress-tested against brute force
 */
#pragma once

template<class T>
vector<T> maxPlusCC(const vector<T>& a, const vector<T>& b) {
	int n = sz(a), m = sz(b), i = 0, j = 0;
	vector<T> c(n + m - 1);
	c[0] = a[0] + b[0];
	while (i + j < n + m - 2) {
		if (j == m-1 || (i < n-1 && a[i+1] - a[i] >= b[j+1] - b[j])) i++;
		else j++;
		c[i + j] = a[i] + b[j];
	}
	return c;
}

template<class T>
vector<T> maxPlusCA(const vector<T>& a, const vector<T>& b) {
	int n = sz(a), m = sz(b);
	vector<T> c(n + m - 1);
	// best j for k in [kl, kr) lies in [jl, jr]
	auto rec = [&](auto& self, int kl, int kr, int jl, int jr) -> void {
		if (kl >= kr) return;
		int k = (kl + kr) / 2, lo = max(jl, k-n+1), hi = min(jr, k), o = lo;
		rep(j,lo,hi+1) if (b[j] + a[k-j] > b[o] + a[k-o]) o = j;
		c[k] = b[o] + a[k-o];
		self(self, kl, k, jl, o);
		self(self, k+1, kr, o, jr);
	};
	rec(rec, 0, n + m - 1, 0, m - 1);
	return c;
}
