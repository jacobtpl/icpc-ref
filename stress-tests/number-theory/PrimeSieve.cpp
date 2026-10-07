#include "../utilities/template.h"

#include "../../content/number-theory/PrimeSieve.h"

// independent oracle: plain Eratosthenes with smallest prime factors
vi spf;
void check() {
	assert(sz(cp) == LIM && sz(nx) == LIM && sz(lp) == LIM && sz(cnt) == LIM);
	vi primes;
	rep(i,2,LIM) if (spf[i] == i) primes.push_back(i);
	assert(pr == primes);
	rep(i,0,2) assert(!cp[i] && lp[i] == -1 && nx[i] == -1 && cnt[i] == -1);
	rep(i,2,LIM) {
		int p = spf[i], m = i, c = 0;
		while (m % p == 0) m /= p, c++;
		assert(cp[i] == (p != i));
		assert(lp[i] >= 0 && lp[i] < sz(pr) && pr[lp[i]] == p);
		assert(nx[i] == m && cnt[i] == c);
	}
}

int main() {
	spf.assign(LIM, 0);
	for (int i = 2; i < LIM; i++) if (!spf[i])
		for (int j = i; j < LIM; j += i) if (!spf[j]) spf[j] = i;
	// brute-force sanity check of the oracle itself
	rep(i,2,20000) {
		int p = 2;
		while (i % p) p++;
		assert(spf[i] == p);
	}
	sieve();
	check();
	// factorize through lp/nx/cnt
	rep(i,2,LIM) if (i % 97 == 0 || i < 100000) {
		ll prod = 1; int last = -1;
		for (int x = i; x > 1; x = nx[x]) {
			assert(lp[x] > last); last = lp[x];
			rep(j,0,cnt[x]) prod *= pr[lp[x]];
		}
		assert(prod == i);
	}
	// sieve() re-initialises all its tables, so calling it again must be harmless
	sieve();
	check();
	cout<<"Tests passed!"<<endl;
}
