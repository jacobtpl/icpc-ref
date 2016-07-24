/**
 * Author: Lukas Polacek
 * Date: 2010-01-28
 * Source: Wikipedia
 * Description: Pollard's rho algorithm. It is a probabilistic factorisation
 * algorithm, whose expected time complexity is good. Before you start using it,
 * run {\tt init(bits)}, where bits is the length of the numbers you use.
 * Time: Expected running time should be good enough for 50-bit numbers.
 */
#pragma once

#include "miller_rabin.h"
#include "eratosthenes.h"
#include "euclid.h"

vector<ull> pr;
ull f(ull a, ull n, ull &has) {
	return (mod_mul(a, a, n) + has) % n;
}
vector<ull> factor(ull d) {
	vector<ull> res;
	for (size_t i = 0; i < pr.size() && pr[i]*pr[i] <= d; i++)
		if (d % pr[i] == 0) {
			while (d % pr[i] == 0) d /= pr[i];
			res.push_back(pr[i]);
		}
	//d is now a product of at most 2 primes.
	if (d > 1) {
		if (prime(d))
			res.push_back(d);
		else while (true) {
			ull has = rand() % 2321 + 47;
			ull x = 2, y = 2, c = 1;
			for (; c==1; c = gcd((y > x ? y - x : x - y), d)) {
				x = f(x, d, has);
				y = f(f(y, d, has), d, has);
			}
			if (c != d) {
				res.push_back(c); d /= c;
				if (d != c) res.push_back(d);
				break;
			}
		}
	}
	return res;
}
void init(int bits) {//how many bits do we use?
	vi p = eratosthenes_sieve(1 << ((bits + 2) / 3));
	vector<ull> pr(p.size());
	for (size_t i=0; i<pr.size(); i++)
		pr[i] = p[i];
}
