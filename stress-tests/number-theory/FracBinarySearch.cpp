#include "../utilities/template.h"

#include "../../content/number-theory/FracBinarySearch.h"

#include "../../content/number-theory/euclid.h"

typedef __int128 LL;
mt19937_64 rng(5);
ll calls, NN;

void testExact() {
	// tiny N, exhaustive oracle, thresholds a/b with both >= and >
	rep(N,1,60) rep(it,0,10000) {
		ll b = (ll)(rng() % (3*N)) + 1, a = (ll)(rng() % (b+1));
		bool strict = rng() % 2;
		if (strict && a == b) continue;
		auto fn = [&](Frac f) {
			assert(f.q >= 1 && f.q <= N && f.p >= 0 && f.p <= N);
			return strict ? f.p*b > a*f.q : f.p*b >= a*f.q;
		};
		Frac r = fracBS(fn, N);
		ll bp = 1, bq = 1;
		for (ll q = 1; q <= N; q++) for (ll p = 0; p <= q; p++)
			if (fn(Frac{p, q}) && p*bq < bp*q) bp = p, bq = q;
		assert(r.p == bp && r.q == bq);
	}
	// large N: threshold p0/q0 must be found exactly, and "> p0/q0" must
	// give its successor in the Farey sequence F_N (computed with euclid).
	// f must never be called with p or q > N, and only O(log N) times.
	for (ll N : {1000LL, 1000000LL, 1000000000LL, 1000000000000LL,
			1000000000000000000LL, 3000000000000000000LL
#ifdef TEST_HUGE_N
			// signed overflow of lo.q * adv for N > 2^63 / 3
			, (1LL << 62) - 1
#endif
			}) {
		rep(it,0,50000) {
			ll q0 = (ll)(rng() % N) + 1, p0 = (ll)(rng() % (q0+1));
			if (it % 7 == 0) q0 = N - (ll)(rng() % min(N, 100LL)), p0 = (ll)(rng() % (q0+1));
			if (it % 11 == 0) p0 = min(q0, (ll)(rng() % 5));
			if (it % 13 == 0) p0 = q0 - min(q0, (ll)(rng() % 5));
			ll g = __gcd(p0, q0); p0 /= g; q0 /= g;
			calls = 0; NN = N;
			Frac r = fracBS([&](Frac f) {
				calls++; assert(f.q >= 1 && f.q <= NN && f.p >= 0 && f.p <= NN);
				return (LL)f.p*q0 >= (LL)p0*f.q; }, N);
			assert(r.p == p0 && r.q == q0);
			assert(calls <= 4 * (64 - __builtin_clzll(N)) + 10);
			if (p0 == q0) continue;
			ll x, y; euclid(p0, q0, x, y);
			ll q1 = ((-x) % q0 + q0) % q0;
			q1 += (N - q1) / q0 * q0;
			ll p1 = (ll)(((LL)p0*q1 + 1) / q0);
			calls = 0;
			r = fracBS([&](Frac f) {
				calls++; assert(f.q >= 1 && f.q <= NN && f.p >= 0 && f.p <= NN);
				return (LL)f.p*q0 > (LL)p0*f.q; }, N);
			assert(r.p == p1 && r.q == q1);
			assert(calls <= 4 * (64 - __builtin_clzll(N)) + 10);
		}
	}
	// usage example from the header: smallest p/q >= 1/3 is {1,3}
	Frac r = fracBS([](Frac f) { return 3*f.p>=f.q; }, 10);
	assert(r.p == 1 && r.q == 3);
#ifdef TEST_OLD_USAGE
	// the example as it used to be written, "f.p>=3*f.q", is false on all
	// of [0, 1] and dies on assert(f(hi)) instead of returning {1,3}
	r = fracBS([](Frac f) { return f.p>=3*f.q; }, 10);
	assert(r.p == 1 && r.q == 3);
#endif
	r = fracBS([](Frac) { return true; }, 10);
	assert(r.p == 0 && r.q == 1);
}

int main() {
	testExact();
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
			double x = rand() / (RAND_MAX + 1.0);
			// x *= min(n, 10); // if testing with search range (0, n)
			auto fn = [&](Frac f) { return (double)f.p >= x * (double)f.q; };
			Frac f = fracBS(fn, n);
			auto best = actual.lower_bound(x)->second;
			assert(best.first == f.p);
			assert(best.second == f.q);
		}
	}
	cout<<"Tests passed!"<<endl;
	return 0;
}
