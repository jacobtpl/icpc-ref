#include "../utilities/template.h"

#include "../../content/number-theory/ModPow.h"

ll oracle(ll b, ll e) {
	__int128 r = 1, x = b % mod;
	for (; e; e /= 2, x = x * x % mod)
		if (e & 1) r = r * x % mod;
	return (ll)r;
}

void check(ll b, ll e) {
	ll r = modpow(b, e), w = oracle(b, e);
	if (r != w) {
		cout << "modpow(" << b << ", " << e << ") = " << r << ", expected " << w << endl;
		exit(1);
	}
}

int main() {
	assert(mod == 1000000007);
	// naive repeated multiplication for small values
	rep(b,0,200) {
		ll r = 1;
		rep(e,0,200) {
			assert(modpow(b, e) == r);
			assert(modpow(mod - b - 1, e) == oracle(mod - b - 1, e));
			r = r * b % mod;
		}
	}
	assert(modpow(0, 0) == 1 && modpow(0, 5) == 0 && modpow(1, LLONG_MAX) == 1);
	// Fermat
	rep(b,1,1000) assert(modpow(b, mod - 1) == 1 && modpow(b, mod - 2) * b % mod == 1);
	mt19937_64 rng(3);
	// b * b must fit in a signed 64-bit integer: 0 <= b <= 3037000499
	const ll BLIM = 3037000499LL;
	rep(it,0,2000000) {
		ll b = it % 4 == 0 ? (ll)(rng() % mod) : it % 4 == 1 ? mod - 1 - (ll)(rng() % 100)
			: it % 4 == 2 ? BLIM - (ll)(rng() % 100) : (ll)(rng() % (BLIM + 1));
		ll e = it % 3 == 0 ? (ll)(rng() % 100) : it % 3 == 1 ? LLONG_MAX - (ll)(rng() % 100)
			: (ll)(rng() >> (1 + rng() % 63));
		check(b, e);
	}
#ifdef BIG_B
	// Opt-in: bases that are not reduced below ~3.04e9 overflow in b * b.
	check(3037000500LL, 2);
	check(10000000000LL, 2);
	check(-2, 3);
#endif
	cout<<"Tests passed!"<<endl;
}
