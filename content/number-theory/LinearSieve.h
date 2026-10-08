/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: folklore (Euler's sieve), https://cp-algorithms.com/algebra/prime-sieve-linear.html
 * Description: Linear sieve. For $2 \le i <$ LIM computes the smallest
 * prime factor $lp[i]$, M\"obius $mu[i]$ and Euler's $phi[i]$; pr holds
 * all primes $<$ LIM in increasing order. Needs LIM $\ge 2$. $lp[0]=lp[1]=0$,
 * $mu[1]=phi[1]=1$, $mu[0]=phi[0]=0$. Every composite $ip$ is visited
 * exactly once, with $p=lp[ip]$, so any multiplicative $f$ fits
 * the same three cases: $f[1]=1$; $f[p]$ for primes; $f[ip]=f[i]f[p]$
 * if $p \nmid i$; and a rule for $p \mid i$. If that rule needs the
 * exponent, also keep $q[i]$ = the largest power of $lp[i]$ dividing $i$
 * ($q[p]=p$, $q[ip]=p$ if $p \nmid i$, else $q[i]p$): then
 * $f[ip]=f[i/q[i]] \cdot f(q[i]p)$.
 * Usage: sieve(); // once
 *  while (x > 1) { int p = lp[x]; x /= p; } // factorise x < LIM
 * Time: O(LIM), LIM=1e7 $\approx$ 0.1s (120 MB), 1e8 $\approx$ 1.3s
 * Status: stress-tested against trial division and naive mu/phi
 */
#pragma once

const int LIM = 1e7;
int lp[LIM], mu[LIM], phi[LIM];
vi pr;

void sieve() {
	mu[1] = phi[1] = 1;
	rep(i,2,LIM) {
		if (!lp[i]) // i is prime
			lp[i] = i, mu[i] = -1, phi[i] = i - 1, pr.pb(i);
		for (int p : pr) {
			if (p > lp[i] || p > (LIM - 1) / i) break;
			int k = i * p;
			lp[k] = p;
			if (p == lp[i]) // p | i
				mu[k] = 0, phi[k] = phi[i] * p;
			else // p, i coprime
				mu[k] = -mu[i], phi[k] = phi[i] * (p - 1);
		}
	}
}
