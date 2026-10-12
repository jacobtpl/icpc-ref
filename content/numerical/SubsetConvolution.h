/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: Bj\"orklund, Husfeldt, Kaski, Koivisto, "Fourier meets M\"obius: fast subset convolution" (2007)
 * Description: \texttt{zeta(a)} replaces $a[S]$ by
 * $\sum_{U \subseteq S} a[U]$; \texttt{zeta(a, 1)} is the inverse
 * (M\"obius transform). \texttt{subsetConv} returns
 * $c[S] = \sum_{U \subseteq S} a[U] \cdot b[S \setminus U]$.
 * Masks are bitsets over bits $0..n-1$; $a$ and $b$ must have the
 * same size $N = 2^n \ge 1$. The value type is a modular type or
 * \texttt{unsigned long long} (exact mod $2^{64}$; intermediate
 * values overflow \texttt{ll} even when the answer fits).
 * \texttt{subsetConv} stores $2N(n+1)$ values; never pass \texttt{w}.
 * Time: O(N \log N) for zeta, O(N \log^2 N) for subsetConv
 * ($n=20$: $\approx 1.2$s with mint, $\approx 1.7$s and 340 MB for 64-bit).
 * Usage: vector<mint> c = subsetConv(a, b);
 * Status: stress-tested
 */
#pragma once

template<class T>
void zeta(vector<T>& a, bool inv = 0, int w = 1) {
	for (int n = sz(a) / w, s = 1; s < n; s *= 2)
		rep(i,0,n) if (i & s) rep(k,0,w) {
			T &x = a[i*w + k], y = a[(i^s)*w + k];
			if (inv) x -= y; else x += y;
		}
}
template<class T>
vector<T> subsetConv(const vector<T>& a, const vector<T>& b) {
	int N = sz(a), w = __builtin_ctz(N) + 1;
	vector<T> f(N * w), g(f), c(N);
	rep(i,0,N) {
		int p = i*w + __builtin_popcount(i);
		f[p] = a[i]; g[p] = b[i];
	}
	zeta(f, 0, w); zeta(g, 0, w);
	rep(i,0,N) for (int k = w; k--;) {
		T x = 0;
		rep(j,0,k+1) x += f[i*w + j] * g[i*w + k - j];
		f[i*w + k] = x;
	}
	zeta(f, 1, w);
	rep(i,0,N) c[i] = f[i*w + __builtin_popcount(i)];
	return c;
}
