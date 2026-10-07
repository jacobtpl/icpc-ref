#include "../utilities/template.h"

#include "../../content/number-theory/SiyongModular.h"
#ifdef TEST_EUCLID_CONFLICT
// Opt-in: SiyongModular.h carries its own copy of euclid(), so pasting it
// together with euclid.h (needed by CRT.h) is a redefinition error.
#include "../../content/number-theory/euclid.h"
#endif

typedef __int128 i128;
typedef unsigned long long ull;
mt19937_64 rng(998244353);

ll red(i128 z) { return (ll)((z % MOD + MOD) % MOD); }
ll mpow(ll a, ull e) {
	ll r = 1;
	for (a = red(a); e; e >>= 1, a = a * a % MOD) if (e & 1) r = r * a % MOD;
	return r;
}
ll rndVal() {
	switch (rng() % 6) {
		case 0: return (ll)(rng() % 5);
		case 1: return MOD - 1 - (ll)(rng() % 5);
		case 2: return MOD + (ll)(rng() % 5) - 2;
		case 3: return (ll)(rng() % MOD);
		case 4: return (ll)rng();
		default: return rng() % 2 ? LLONG_MAX - (ll)(rng() % 3) : LLONG_MIN + (ll)(rng() % 3);
	}
}
template<class T> void ctor(T z) {
	mint a(z);
	assert(0 <= a.v && a.v < MOD);
	assert(a.v == red((i128)z));
}

int main() {
	assert(mint().v == 0);
	// construction from every integer type, including negative / out of range values
	for (i128 z : {(i128)0, (i128)1, (i128)-1, (i128)MOD, (i128)MOD - 1, (i128)MOD + 1,
			-(i128)MOD, -(i128)MOD - 1, -(i128)MOD + 1, (i128)2 * MOD, (i128)INT_MAX, (i128)INT_MIN,
			(i128)UINT_MAX, (i128)LLONG_MAX, (i128)LLONG_MIN, (i128)ULLONG_MAX, (i128)SHRT_MIN,
			(i128)SHRT_MAX, (i128)-128, (i128)127, (i128)255, (i128)65535, (i128)-2, (i128)-12345}) {
		ctor(z); ctor(-z); ctor((i128)(ll)z * (ll)z);
		ctor((int)z); ctor((unsigned)z); ctor((ll)z); ctor((ull)z); ctor((long)z);
		ctor((unsigned __int128)z >> 1);
		ctor((short)z); ctor((unsigned short)z);
		ctor((signed char)z); ctor((unsigned char)z); ctor((char)z); ctor((bool)z);
	}
	rep(it,0,200000) {
		ll z = rndVal();
		ctor(z); ctor((int)z); ctor((unsigned)z); ctor((ull)z); ctor((short)z); ctor((signed char)z);
		ctor((i128)z * rndVal());
	}
	// arithmetic against a 128-bit oracle
	rep(it,0,1000000) {
		ll x = rndVal(), y = rndVal();
		mint a(x), b(y);
		assert((a + b).v == red((i128)x + y));
		assert((a - b).v == red((i128)x - y));
		assert((a * b).v == red((i128)red(x) * red(y)));
		assert((-a).v == red(-(i128)x));
		assert((int)a == a.v);
		mint c = a; c += b; assert(c.v == (a + b).v);
		c = a; c -= b; assert(c.v == (a - b).v);
		c = a; c *= b; assert(c.v == (a * b).v);
		// mixed with plain integers on either side
		assert((a + y).v == (a + b).v && (x + b).v == (a + b).v);
		assert((a - y).v == (a - b).v && (x - b).v == (a - b).v);
		assert((a * y).v == (a * b).v && (x * b).v == (a * b).v);
		// aliasing
		c = a; c += c; assert(c.v == red(2 * (i128)red(x)));
		c = a; c *= c; assert(c.v == red((i128)red(x) * red(x)));
		c = a; c -= c; assert(c.v == 0);
		if (b.v) {
			mint i = invert(b);
			assert(0 <= i.v && i.v < MOD && (i * b).v == 1);
			assert(i.v == mpow(y, MOD - 2));
			mint q = a / b;
			assert((q * b).v == a.v && q.v == red((i128)red(x) * i.v));
			c = a; c /= b; assert(c.v == q.v);
			c = b; c /= c; assert(c.v == 1);
			assert((x / b).v == q.v && (a / y).v == q.v);
		}
		if (it % 8 == 0) {
			ull e = it % 16 ? rng() : rng() % 40;
			assert(pow(a, e).v == mpow(x, e));
			assert(pow(a, (ll)(e >> 1)).v == mpow(x, e >> 1));
			assert(pow(a, (int)(e >> 34)).v == mpow(x, e >> 34));
		}
	}
	// const objects
	const mint k(5);
	assert((-k).v == MOD - 5 && (int)k == 5 && (k + k).v == 10 && pow(k, 2).v == 25);
	assert(pow(mint(0), 0).v == 1 && pow(mint(0), 5).v == 0 && pow(mint(7), 1).v == 7);
	assert(pow(mint(3), MOD - 1).v == 1 && pow(mint(3), (MOD - 1) / 2).v == MOD - 1);
	// exhaustive-ish inverses of small and large residues
	rep(i,1,200000) {
		assert((invert(mint(i)) * i).v == 1);
		assert((invert(mint(MOD - i)) * (MOD - i)).v == 1);
	}
	// long chains stay normalised
	mint s, f(1);
	ll rs = 0, rf = 1;
	rep(i,1,2000000) {
		s += f; f *= i; s -= mint(i) * i;
		rs = red((i128)rs + rf); rf = rf * i % MOD; rs = red(rs - (ll)i * i);
		assert(0 <= s.v && s.v < MOD);
	}
	assert(s.v == rs && f.v == rf);
	cout<<"Tests passed!"<<endl;
}
