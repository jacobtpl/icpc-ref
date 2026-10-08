/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: Bj\"orklund, Husfeldt, Kaski, Koivisto, "Fourier meets M\"obius: fast subset convolution" (2007)
 * Description: \texttt{zeta(a)} replaces $a[S]$ by
 * $\sum_{T \subseteq S} a[T]$; \texttt{zeta(a, 1)} is the inverse
 * (M\"obius transform). \texttt{subsetConv} returns
 * $c[S] = \sum_{T \subseteq S} a[T] \cdot b[S \setminus T]$.
 * Masks are bitsets over bits $0..n-1$; $a$ and $b$ must have the
 * same size $N = 2^n \ge 1$. $T$ is a modular type or
 * \texttt{unsigned ll}: intermediate values overflow \texttt{ll}
 * even when the answer fits, but the result is exact mod $2^{64}$.
 * \texttt{subsetConv} stores $2N(n+1)$ values
 * ($n=20$: 170 MB for mint, 340 MB for 64-bit).
 * Leave \texttt{w} (internal block width) as 1.
 * Time: O(N \log N) for zeta, O(N \log^2 N) for subsetConv
 * ($n=20$: $\approx 1$s with mint).
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
