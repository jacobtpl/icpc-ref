#include "../utilities/template.h"

#include "../../content/numerical/ShittyBigUInt.h"

typedef unsigned __int128 u128;
mt19937_64 rng(12345);

template<int B> u128 val(const BigUInt<B>& a) {
	u128 r = 0;
	for (int i = sz(a); i--;) r = (r << B) | (u128)a[i];
	return r;
}
template<int B> BigUInt<B> make(u128 v) {
	BigUInt<B> r;
	for (; v; v >>= B) r.push_back((int)(v & ((1 << B) - 1)));
	return r;
}
template<int B> void canon(const BigUInt<B>& a) {
	assert(a.empty() || a.back() > 0);
	for (int x : a) assert(0 <= x && x < (1 << B));
}
// random value with a random bit length, biased towards limbs that are all-ones / zero
template<int B> u128 rnd(int maxBits) {
	int bits = (int)(rng() % (maxBits + 1));
	u128 v = ((u128)rng() << 64) | rng();
	int mode = (int)(rng() % 4);
	if (mode == 0) v = ~(u128)0;
	if (mode == 1) v = (u128)1 << (rng() % 127), v -= rng() % 2;
	if (mode == 2) for (int i = 0; i * B + B <= 128; i++) if (rng() % 2) v |= (u128)((1 << B) - 1) << (i * B);
	if (bits == 0) return 0;
	if (bits < 128) v &= (((u128)1 << bits) - 1);
	return v;
}

template<int B> void testSmall(int iters) {
	typedef BigUInt<B> Big;
	rep(it,0,iters) {
		u128 x = rnd<B>(100), y = rnd<B>(100);
		Big a = make<B>(x), b = make<B>(y);
		assert(val(a) == x && val(b) == y);
		{ Big c = a; c += b; canon(c); assert(val(c) == x + y); }
		{ Big c = a + b; canon(c); assert(val(c) == x + y); }
		{ Big c = a; c += c; canon(c); assert(val(c) == x + x); }
		if (x >= y) {
			Big c = a; c -= b; canon(c); assert(val(c) == x - y);
			Big d = a - b; canon(d); assert(val(d) == x - y);
		}
		assert((a < b) == (x < y));
		assert((a > b) == (x > y));
		unsigned s = (unsigned)(rng() % 27);
		{ Big c = a << s; canon(c); assert(val(c) == x << s); }
		{ Big c = a; c <<= s; canon(c); assert(val(c) == x << s); }
		s = (unsigned)(rng() % 110);
		{ Big c = a >> s; canon(c); assert(val(c) == x >> s); }
		{ Big c = a; c >>= s; canon(c); assert(val(c) == x >> s); }
	}
	// exhaustive tiny values
	rep(x,0,70) rep(y,0,70) {
		Big a(x), b(y);
		canon(a); assert(val(a) == (u128)x);
		assert(val(a + b) == (u128)(x + y));
		if (x >= y) assert(val(a - b) == (u128)(x - y));
		assert((a < b) == (x < y) && (a > b) == (x > y));
		rep(s,0,12) {
			assert(val(a << (unsigned)s) == (u128)x << s);
			assert(val(a >> (unsigned)s) == (u128)x >> s);
		}
	}
	// ll constructor
	rep(it,0,2000) {
		ll x = (ll)(rng() >> (1 + rng() % 63));
		Big a(x); canon(a); assert(val(a) == (u128)x);
	}
	{ Big a(LLONG_MAX); canon(a); assert(val(a) == (u128)LLONG_MAX); }
}

// long numbers: algebraic identities
template<int B> void testLong(int iters, int maxLimbs) {
	typedef BigUInt<B> Big;
	auto gen = [&](int n) {
		Big a;
		int mode = (int)(rng() % 3);
		rep(i,0,n) a.push_back(mode == 0 ? (1 << B) - 1 : mode == 1 ? (int)(rng() % 2) * ((1 << B) - 1) : (int)(rng() % (1 << B)));
		while (!a.empty() && a.back() == 0) a.pop_back();
		return a;
	};
	rep(it,0,iters) {
		Big a = gen((int)(rng() % (maxLimbs + 1))), b = gen((int)(rng() % (maxLimbs + 1)));
		Big s = a + b, s2 = b + a;
		canon(s); assert(s == s2);
		assert(!(s < a) && !(s < b));
		Big d = s - b; canon(d); assert(d == a);
		d = s - a; canon(d); assert(d == b);
		unsigned sh = (unsigned)(rng() % (3 * B + 1));
		Big l = a << sh; canon(l);
		Big r = l >> sh; canon(r); assert(r == a);
		// (a >> sh) << sh + low bits == a
		Big hi = a >> sh; canon(hi);
		Big back = hi << sh; canon(back);
		assert(!(a < back));
		Big low = a - back; canon(low);
		Big lim = Big(1) << sh;
		assert(low < lim);
		// doubling equals shift by one
		assert(a + a == a << 1);
	}
}

int main() {
	testSmall<20>(100000);
	testSmall<15>(50000);
	testSmall<7>(50000);
	testSmall<1>(5000);
	testLong<20>(20000, 40);
	testLong<7>(20000, 40);
	cout<<"Tests passed!"<<endl;
}
