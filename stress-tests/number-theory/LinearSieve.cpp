#include "../utilities/template.h"
#define pb push_back
#include "../../content/number-theory/LinearSieve.h"

// Copy of the header's sieve with LIM as a parameter, for small / edge limits.
namespace small {
template <int L> struct S {
	vi lp, mu, phi, pr;
	S() : lp(L), mu(L), phi(L) {
		mu[1] = phi[1] = 1;
		rep(i,2,L) {
			if (!lp[i]) lp[i] = i, mu[i] = -1, phi[i] = i - 1, pr.pb(i);
			for (int p : pr) {
				if (p > lp[i] || p > (L - 1) / i) break;
				int k = i * p;
				lp[k] = p;
				if (p == lp[i]) mu[k] = 0, phi[k] = phi[i] * p;
				else mu[k] = -mu[i], phi[k] = phi[i] * (p - 1);
			}
		}
	}
};
}

struct Naive { int lp, mu, phi; };
Naive naive(int n) {
	if (n < 2) return {0, n, n};
	Naive r{0, 1, n};
	int x = n;
	for (int p = 2; (ll)p * p <= x; p++) if (x % p == 0) {
		if (!r.lp) r.lp = p;
		int c = 0;
		while (x % p == 0) x /= p, c++;
		r.mu = c > 1 ? 0 : -r.mu;
		r.phi = r.phi / p * (p - 1);
	}
	if (x > 1) {
		if (!r.lp) r.lp = x;
		r.mu = -r.mu;
		r.phi = r.phi / x * (x - 1);
	}
	return r;
}

void check(int n, int l, int m, int ph) {
	Naive e = naive(n);
	if (l != e.lp || m != e.mu || ph != e.phi) {
		cout << "differ at n=" << n << ": lp " << l << " vs " << e.lp
		     << ", mu " << m << " vs " << e.mu << ", phi " << ph << " vs " << e.phi << endl;
		assert(false);
	}
}

template <int L> void checkSmall() {
	small::S<L> s;
	vi pr;
	rep(n,0,L) {
		check(n, s.lp[n], s.mu[n], s.phi[n]);
		if (n >= 2 && naive(n).lp == n) pr.pb(n);
	}
	assert(s.pr == pr);
}

int main() {
	sieve();

	// primes vs. plain Eratosthenes, lp[i]==i iff prime
	vector<bool> comp(LIM);
	vi pr2;
	rep(i,2,LIM) if (!comp[i]) {
		pr2.pb(i);
		for (ll j = (ll)i * i; j < LIM; j += i) comp[j] = 1;
	}
	assert(pr == pr2);
	assert(lp[0] == 0 && lp[1] == 0 && mu[0] == 0 && phi[0] == 0 && mu[1] == 1 && phi[1] == 1);
	rep(i,2,LIM) {
		assert((lp[i] == i) == !comp[i]);
		assert(lp[i] >= 2 && i % lp[i] == 0 && !comp[lp[i]]);
		assert(mu[i] >= -1 && mu[i] <= 1 && phi[i] >= 1 && phi[i] < i);
	}

	// exhaustive vs. trial division for small n, random for large n
	rep(n,0,300000) check(n, lp[n], mu[n], phi[n]);
	for (int n = LIM - 1; n >= LIM - 1000; n--) check(n, lp[n], mu[n], phi[n]);
	mt19937 rng(1);
	rep(it,0,1000000) {
		int n = (int)(rng() % LIM);
		check(n, lp[n], mu[n], phi[n]);
	}

	// phi by definition
	rep(n,1,2001) {
		int c = 0;
		rep(k,1,n+1) c += __gcd(n, k) == 1;
		assert(phi[n] == c);
	}

	// sum_{d|n} mu(d) = [n==1], sum_{d|n} phi(d) = n
	const int N = 200000;
	vi smu(N), sphi(N);
	rep(d,1,N) for (int n = d; n < N; n += d) smu[n] += mu[d], sphi[n] += phi[d];
	rep(n,1,N) assert(smu[n] == (n == 1) && sphi[n] == n);

	// factorisation from the usage comment
	rep(it,0,200000) {
		int n = it < 1000 ? it + 1 : (int)(rng() % (LIM - 1)) + 1;
		int x = n; ll prod = 1; int last = 1;
		while (x > 1) {
			int p = lp[x]; x /= p;
			assert(p >= last && lp[p] == p);
			last = p; prod *= p;
		}
		assert(prod == n);
	}

	// the q[] recipe from the description, applied to the divisor count
	{
		vi q(N), f(N), cnt(N);
		f[1] = 1;
		rep(i,2,N) {
			int p = lp[i];
			if (p == i) q[i] = p, f[i] = 2;
			else {
				int j = i / p;
				if (lp[j] != p) q[i] = p, f[i] = f[j] * 2;
				else {
					q[i] = q[j] * p;
					int e = 0; for (int t = q[i]; t > 1; t /= p) e++;
					f[i] = f[j / q[j]] * (e + 1);
				}
			}
		}
		rep(d,1,N) for (int n = d; n < N; n += d) cnt[n]++;
		rep(n,1,N) assert(f[n] == cnt[n]);
	}

	// small and awkward limits, incl. the i*p overflow boundary
	checkSmall<2>(); checkSmall<3>(); checkSmall<4>(); checkSmall<5>();
	checkSmall<6>(); checkSmall<7>(); checkSmall<8>(); checkSmall<9>();
	checkSmall<10>(); checkSmall<11>(); checkSmall<12>(); checkSmall<13>();
	checkSmall<16>(); checkSmall<17>(); checkSmall<25>(); checkSmall<27>();
	checkSmall<30>(); checkSmall<32>(); checkSmall<49>(); checkSmall<50>();
	checkSmall<64>(); checkSmall<100>(); checkSmall<121>(); checkSmall<122>();
	checkSmall<1000>(); checkSmall<46340>(); checkSmall<46341>();
	checkSmall<46342>(); checkSmall<65536>(); checkSmall<65537>();

	cout << "Tests passed!" << endl;
}
