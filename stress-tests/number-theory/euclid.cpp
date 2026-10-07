#include "../utilities/template.h"

#include "../../content/number-theory/euclid.h"

typedef __int128 i128;
mt19937_64 rng(2024);
ll absll(ll v) { return v < 0 ? -v : v; }

// checks a*x + b*y = g, |g| = gcd(|a|,|b|) and that the coefficients are small
void check(ll a, ll b, bool nonNeg) {
	ll x = 1234567, y = 7654321, g = euclid(a, b, x, y);
	ll G = __gcd(absll(a), absll(b));
	assert((i128)a * x + (i128)b * y == g);
	assert(absll(g) == G);
	if (nonNeg) assert(g == G);
	if (G) {
		assert(absll(x) <= max(1LL, absll(b) / G));
		assert(absll(y) <= max(1LL, absll(a) / G));
	}
	if (nonNeg && G == 1 && b > 0) // x is the inverse of a mod b
		assert((ll)((((i128)a * x) % b + b) % b) == 1 % b);
}

ll rndBits(int bits) {
	int b = (int)(rng() % (bits + 1));
	return b ? (ll)(rng() & ((1ULL << b) - 1)) : 0;
}

int main() {
	rep(a,-150,151) rep(b,-150,151) check(a, b, a >= 0 && b >= 0);
	rep(it,0,2000000) {
		ll a = rndBits(62), b = rndBits(62);
		check(a, b, true);
		if (it % 4 == 0) {
			ll g = rndBits(30) + 1; // force a large common factor
			check(a / g * g, b / g * g, true);
		}
		if (rng() % 2) a = -a;
		if (rng() % 2) b = -b;
		check(a, b, false);
	}
	// sign convention: the result is gcd(a,b) >= 0 for a,b >= 0, but can be -gcd otherwise
	{ ll x, y; assert(euclid(-4, 6, x, y) == 2 && euclid(4, -6, x, y) == -2 && euclid(0, -5, x, y) == -5); }
	// extremes
	const ll M = LLONG_MAX;
	for (ll a : {0LL, 1LL, 2LL, M, M - 1, M / 2, M / 2 + 1, (ll)1e18, (1LL << 62)})
		for (ll b : {0LL, 1LL, 2LL, M, M - 1, M / 2, M / 2 + 1, (ll)1e18, (1LL << 62)})
			check(a, b, true), check(-a, b, false), check(a, -b, false), check(-a, -b, false);
	// worst case recursion depth: consecutive Fibonacci numbers
	ll f0 = 0, f1 = 1;
	while (f1 <= M - f0) {
		ll t = f0 + f1; f0 = f1; f1 = t;
		check(f1, f0, true); check(f0, f1, true);
	}
	cout<<"Tests passed!"<<endl;
}
