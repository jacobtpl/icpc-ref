/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: own work; idea from Claude, Navarro, Ordonez, "The wavelet matrix" (2015)
 * Description: Wavelet matrix on a static array $a$ (any ordered type, values are compressed).
 *  kth(l, r, k) returns the $k$-th smallest ($0$-indexed) element of $a[l, r)$, requires $0 \le k < r - l$.
 *  count(l, r, v) returns the number of $i$ in $[l, r)$ with $a[i] < v$, requires $0 \le l \le r \le N$.
 *  The number of values in $[x, y)$ is count(l, r, y) - count(l, r, x).
 * Usage:
 *  Wavelet<int> w(a);
 *  w.kth(l, r, (r - l) / 2); // upper median of a[l, r)
 * Time: $O(N \log N)$ build, $O(\log N)$ per query. Memory: $O(N \log N)$ ints.
 * Status: stress-tested
 */
#pragma once

template<class T>
struct Wavelet {
	vector<T> s; vector<vi> b; // b[h][i] = #zeros of bit h in first i
	Wavelet(const vector<T>& a) : s(a) {
		sort(all(s)); s.erase(unique(all(s)), s.end());
		int n = sz(a), L = 0;
		while (sz(s) >> L) L++;
		vi c(n);
		rep(i,0,n) c[i] = (int)(lower_bound(all(s), a[i]) - s.begin());
		b.assign(L, vi(n + 1));
		for (int h = L; h--;) {
			rep(i,0,n) b[h][i+1] = b[h][i] + !(c[i] >> h & 1);
			stable_partition(all(c), [&](int x) { return !(x >> h & 1); });
		}
	}
	T kth(int l, int r, int k) {
		int x = 0;
		for (int h = sz(b); h--;) {
			int p = b[h][l], q = b[h][r], z = b[h].back();
			if (k < q - p) l = p, r = q;
			else k -= q - p, x |= 1 << h, l += z - p, r += z - q;
		}
		return s[x];
	}
	int count(int l, int r, T v) {
		int x = (int)(lower_bound(all(s), v) - s.begin()), res = 0;
		for (int h = sz(b); h--;) {
			int p = b[h][l], q = b[h][r], z = b[h].back();
			if (x >> h & 1) res += q - p, l += z - p, r += z - q;
			else l = p, r = q;
		}
		return res;
	}
};
