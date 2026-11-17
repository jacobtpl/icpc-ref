/**
 * Author: Håkan Terelius
 * Date: 2009-09-25
 * License: CC0
 * Source: http://en.wikipedia.org/wiki/Euler's_totient_function
 * Description: Euler's phi function is defined as phi(n):=# of positive integers <= n that are coprime with n.
 * phi(1)=1, p prime => phi(p^k)=(p-1)p^(k-1), m,n coprime => phi(mn)=phi(m)phi(n).
 * If n=p_1^k_1p_2^k_2 ... p_r^k_r then phi(n) = (p_1-1)p_1^(k_1-1)...(p_r-1)p_r^(k_r-1).
 * phi(n)=n * prod_(p|n)(1-1/p).
 *
 * sum_(d|n) phi(d) = n, sum_(1<= k <= n, gcd(k,n)=1) k = n phi(n)/2, n>1
 *
 * Euler's thm: a,n coprime => a^(phi(n)) == 1 (mod n).
 *
 * Fermat's little thm: p prime => a^(p-1) == 1 (mod p) for all a.
 * Status: Tested
 */
#pragma once

const int LIM = 5000000;
int phi[LIM];

void calculatePhi() {
	rep(i,0,LIM) phi[i] = i&1 ? i : i/2;
	for (int i = 3; i < LIM; i += 2) if(phi[i] == i)
		for (int j = i; j < LIM; j += i) phi[j] -= phi[j] / i;
}
