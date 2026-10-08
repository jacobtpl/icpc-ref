/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: Lucy\_Hedgehog's post in the Project Euler problem 10 thread; own implementation
 * Description: Counts the primes (sums them if $k=1$) up to every distinct
 * value of $\lfloor n/i \rfloor$, for $1 \le n \le 10^{13}$. With
 * $s = \lfloor\sqrt n\rfloor$, \texttt{lo[x]} is the answer for $x \le s$ and
 * \texttt{hi[i]} the answer for $\lfloor n/i \rfloor$, $1 \le i \le s$, so
 * $\pi(n)$ = \texttt{hi[1]}. Counts fit in ll. For sums with $n > 4 \cdot 10^9$
 * use \texttt{T = \_\_int128} (or an unsigned/mod type for the sum modulo).
 * Time: $O(n^{3/4})$, memory $O(\sqrt n)$. $n=10^{11}$ in 0.15s, $n=10^{12}$ in 0.7s.
 * Usage: PrimeCount pc(n); pc.hi[1]; pc.get(n / 7);
 * Status: stress-tested
 */
#pragma once

typedef ll T;
struct PrimeCount {
	ll n; int s; vector<T> lo, hi;
	T f(ll x, bool k) { return k ? T((__int128)x*(x+1)/2-1) : T(x-1); }
	PrimeCount(ll n, bool k = 0) : n(n), s((int)sqrtl(n)), lo(s+1), hi(s+1) {
		rep(i,1,s+1) lo[i] = f(i, k), hi[i] = f(n / i, k);
		for (ll p = 2; p <= s; p++) if (lo[p] != lo[p-1]) {
			T c = lo[p-1], w = k ? p : 1;
			ll m = n / p;
			int a = int(s / p), e = (int)min((ll)s, m / p);
			rep(i,1,a+1) hi[i] -= (hi[i*p] - c) * w;
			rep(i,a+1,e+1) hi[i] -= (lo[ll((double)m/i)] - c) * w;
			for (int i = s; i >= p*p; i--) lo[i] -= (lo[i/p] - c) * w;
		}
	}
	T get(ll x) { return x <= s ? lo[x] : hi[n/x]; } // x = n/i
};
