#include "../utilities/template.h"

#include "../../content/number-theory/ModSqrt.h"

ll mpow(ll b, ll e, ll m) {
	__int128 r = 1 % m, x = b % m;
	for (; e; e /= 2, x = x * x % m)
		if (e & 1) r = r * x % m;
	return (ll)r;
}

bool isPrimeSlow(ll p) {
	if (p < 2) return false;
	for (ll i = 2; i * i <= p; i++) if (p % i == 0) return false;
	return true;
}

void check(ll a, ll p) {
	ll x = sqrt(a, p);
	ll r = ((a % p) + p) % p;
	if (!(0 <= x && x < p) || (ll)((__int128)x * x % p) != r) {
		cout << "sqrt(" << a << ", " << p << ") = " << x << endl;
		exit(1);
	}
}

void testPrime(ll p, int iters, mt19937_64& rng) {
	check(0, p); check(1, p); check(p, p); check(-p, p); check(p + 1, p); check(1 - p, p);
	rep(it,0,iters) {
		// squares are exactly the valid inputs
		ll y = it < 5 ? p - 1 - it : (ll)(rng() % p), a = (ll)((__int128)y * y % p);
		check(a, p);
		check(a - p, p); // negative representative
		check(a + p * (ll)(rng() % 3), p);
	}
}

int main() {
	// exhaustive: all primes below 5000, all residues
	rep(p,2,5000) {
		if (!isPrimeSlow(p)) continue;
		vector<bool> isSq(p);
		rep(y,0,p) isSq[y * y % p] = true;
		rep(a,0,p) {
			if (!isSq[a]) continue; // precondition: a is a square (asserted)
			ll x = sqrt(a, p);
			assert(0 <= x && x < p);
			assert(x * x % p == a);
			if (a < 50) check(a - p, p), check(a + 5LL * p, p);
		}
	}
	mt19937_64 rng(17);
	// large primes with high 2-adic valuation of p-1 (worst case for
	// Tonelli-Shanks) and of every residue class mod 8
	vector<ll> ps = {998244353, 1000000007, 1000000009, 2013265921, 469762049, 167772161,
		754974721, 1004535809, 2147483647, 2113929217, 2281701377LL, 2483027969LL, 2885681153LL,
		65537, 786433, 7340033, 23068673, 104857601, 3037000493LL, 3037000453LL, 2717908993LL, 2130706433LL};
	for (ll p = 3037000499LL; sz(ps) < 60; p--) if (isPrimeSlow(p)) ps.push_back(p);
	while (sz(ps) < 160) {
		ll p = (ll)(rng() % 3037000499LL) + 2;
		if (isPrimeSlow(p)) ps.push_back(p);
	}
	for (ll p : ps) {
		assert(isPrimeSlow(p) && p <= 3037000499LL);
		testPrime(p, 300, rng);
	}
#ifdef BIG_P
	// Opt-in: t * t % p etc. overflow a signed 64-bit integer once p exceeds
	// about 3.04e9 (3221225473 = 3 * 2^30 + 1 is prime).
	for (ll p : {3221225473LL, 1000000000000000003LL, 4611686018427387847LL, 4179340454199820289LL}) {
		rep(it,0,2000) {
			ll y = (ll)(rng() % p), a = (ll)((__int128)y * y % p);
			check(a, p);
		}
	}
#endif
	cout<<"Tests passed!"<<endl;
}
