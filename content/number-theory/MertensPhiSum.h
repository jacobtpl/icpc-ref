/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: own work; idea from https://codeforces.com/blog/entry/54150
 * Description: Prefix sums of multiplicative functions for large $n$ (Du's sieve).
 * \texttt{M(n)} $= \sum_{i=1}^{n} \mu(i)$ from
 * $M(n) = 1 - \sum_{d=2}^{n} M(\lfloor n/d \rfloor)$, and
 * \texttt{Phi(n)} $= \sum_{i=1}^{n} \phi(i) = \sum_{d=1}^{n} \mu(d) T(\lfloor n/d \rfloor)$,
 * $T(m) = m(m+1)/2$. Both are exact for any $n \ge 0$; \texttt{Phi(n)} only
 * fits in a \texttt{ll} for $n \le 5.5 \cdot 10^9$. Call \texttt{initMu()} once first.
 * In general, if $g$ and $f * g$ have easy prefix sums $G$, $H$ then
 * $g(1) S_f(n) = H(n) - \sum_{d=2}^{n} g(d) S_f(\lfloor n/d \rfloor)$.
 * Time: \texttt{initMu} is $O(\texttt{LIM})$; with $\texttt{LIM} \approx n^{2/3}$,
 * \texttt{M(n)} is $O(n^{2/3})$ and \texttt{Phi(n)} is $O(\sqrt n)$ more.
 * $n=10^{10}$ (\texttt{LIM}=5e6) $\approx$ 0.3s, $n=10^{11}$ (\texttt{LIM}=2e7) $\approx$ 1.5s.
 * Usage: initMu(); ll m = M(n); lll s = Phi(n);
 * Status: stress-tested
 */
#pragma once

typedef __int128 lll;
const int LIM = 5e6; // ~ n^(2/3), 4*LIM bytes
int mu[LIM]; // prefix sums of mu after initMu()
map<ll, ll> memo;

void initMu() {
	vi pr; vector<bool> comp(LIM);
	mu[1] = 1;
	rep(i,2,LIM) {
		if (!comp[i]) pr.push_back(i), mu[i] = -1;
		for (int p : pr) {
			if ((ll)i * p >= LIM) break;
			comp[i * p] = 1;
			if (i % p == 0) break;
			mu[i * p] = -mu[i];
		}
	}
	rep(i,1,LIM) mu[i] += mu[i - 1];
}
ll M(ll n) {
	if (n < LIM) return mu[n];
	if (memo.count(n)) return memo[n];
	ll r = 1;
	for (ll l = 2, h; l <= n; l = h + 1)
		h = n / (n / l), r -= (h - l + 1) * M(n / l);
	return memo[n] = r;
}
lll Phi(ll n) {
	lll r = 0;
	for (ll l = 1, h, p = 0, c; l <= n; l = h + 1, p = c) {
		ll q = n / l; h = n / q; c = M(h);
		r += (lll)q * (q + 1) / 2 * (c - p);
	}
	return r;
}
