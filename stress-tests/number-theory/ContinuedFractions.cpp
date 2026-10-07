#include "../utilities/template.h"

#include "../../content/number-theory/ContinuedFractions.h"

typedef long double ld;
mt19937_64 rng(5);

// O(N) oracle: best error over all p/q with 0 <= p <= N, 1 <= q <= N
ld bestErr(double x, ll N) {
	ld best = 1e30L;
	for (ll q = 1; q <= N; q++) {
		ll p0 = (ll)floorl((ld)x * (ld)q);
		for (ll p = p0; p <= p0 + 1; p++)
			best = min(best, fabsl((ld)min(p, N) / (ld)q - (ld)x));
	}
	return best;
}

void testLarge() {
	for (ll N : {1LL, 2LL, 3LL, 7LL, 50LL, 1000LL, 5000LL, 30000LL})
	rep(it,0,N < 100 ? 20000 : 300) {
		double x;
		int t = int(rng() % 5);
		if (t == 0) x = (double)(rng() % 1000000) / 1e6 * (double)(rng() % 3 + 1);
		else if (t == 1) x = (double)(rng() % (3*N+1)) / (double)(rng() % (3*N) + 1);
		else if (t == 2) x = (double)(rng() % 1000000) / 1e6;
		else if (t == 3) x = (double)(rng() % (2*N+2)) + (double)(rng() % 1000) / 1e3;
		else x = (double)(rng() % (2*N+2)); // integers, including 0 and > N
		auto pa = approximate(x, N);
		assert(0 <= pa.first && pa.first <= N && 1 <= pa.second && pa.second <= N);
		assert(__gcd(pa.first, pa.second) == 1);
		ld e = fabsl((ld)pa.first / (ld)pa.second - (ld)x);
		assert(e <= bestErr(x, N) * (1 + 1e-9L) + 1e-17L);
		// documented error bound
		if (x <= 1) assert(e * (ld)N <= 1 + 1e-12L);
#ifdef TEST_QN_BOUND
		// the header used to claim |p/q - x| <= 1/qN, which is false,
		// e.g. N=3, x=0.21: answer 1/3, error 0.123 > 1/9
		if (x <= 1) assert(e * (ld)N * (ld)pa.second <= 1 + 1e-12L);
#endif
	}
	assert(approximate(0, 10) == make_pair(0LL, 1LL));
	assert(approximate(1e15, 10) == make_pair(10LL, 1LL));
	assert(approximate(0.21, 3) == make_pair(1LL, 3LL));
	// "double for N ~ 1e7": exact rationals with p, q <= N must be recovered
	rep(it,0,1000000) {
		ll N = 10000000, q = (ll)(rng() % N) + 1, p = (ll)(rng() % N) + 1, g = __gcd(p, q);
		p /= g, q /= g;
		assert(approximate((double)p / (double)q, N) == make_pair(p, q));
	}
}

int main() {
	testLarge();
	rep(n,1,300) {
		vector<pair<double, pii>> v;
		rep(i,0,n+1) rep(j,1,n+1) if (__gcd(i,j) == 1) {
			double r = (double)i / j;
			v.emplace_back(r, pii(i,j));
		}
		v.emplace_back(1e9, pii(0,0));
		sort(all(v));
		map<double, pii> actual(all(v));

		rep(iter,0,50000) {
			double x = rand() / (RAND_MAX + 1.0) * 3;
			if (rand() % 2 == 0) x = (rand() % (3*n)) / (double)(rand() % (3*n) + 1);
			auto pa = approximate(x, n);
			auto it = actual.lower_bound(x), it2 = it;
			if (it2 != actual.begin()) --it2;
			auto best =
				min(make_pair(abs(it2->first - x), it2->second),
					make_pair(abs(it->first - x), it->second)).second;
			assert(best.first == pa.first);
			assert(best.second == pa.second);
		}
	}
	cout<<"Tests passed!"<<endl;
	return 0;
}
