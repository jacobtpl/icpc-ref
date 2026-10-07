#include "../utilities/template.h"

const ll mod = 1000000007;

#include "../../content/numerical/FastFourierTransformMod.h"

template<ll M> vl simpleConv(vl a, vl b) {
	if (a.empty() || b.empty()) return {};
	int s = sz(a) + sz(b) - 1;
	vl c(s);
	rep(i,0,sz(a)) rep(j,0,sz(b))
		c[i+j] = (c[i+j] + (ll)a[i] * b[j]) % M;
	for(auto &x: c) if (x < 0) x += M;
	return c;
}

int ra() {
	static unsigned X;
	X *= 123671231;
	X += 1238713;
	X ^= 1237618;
	return (X >> 1);
}

mt19937_64 rng(99);
ll rnd(ll lo, ll hi) { return lo + (ll)(rng() % (unsigned long long)(hi - lo + 1)); }

// Small sizes against the naive product, for modulus M. type 0: random,
// 1: all M-1, 2: values in {0, M-1}, 3: values just around multiples of cut.
template<int M> void testSmall(int iters, int maxn) {
	int cut = int(sqrt(M));
	rep(it,0,iters) {
		vl a(rnd(0, maxn)), b(rnd(0, maxn));
		int type = (int)rnd(0, 3);
		auto val = [&]() -> ll {
			if (type == 0) return rnd(0, M - 1);
			if (type == 1) return M - 1;
			if (type == 2) return rnd(0, 1) * (M - 1);
			ll x = rnd(0, cut) * cut + rnd(-1, 1);
			return min(max(x, 0LL), (ll)M - 1);
		};
		for (auto& x : a) x = val();
		for (auto& x : b) x = val();
		assert(simpleConv<M>(a, b) == convMod<M>(a, b));
	}
}

ll mulM(ll a, ll b) { return a * b % mod; }
ll evalM(const vl& v, ll x) {
	ll r = 0;
	for (int i = sz(v) - 1; i >= 0; i--) r = (mulM(r, x) + v[i]) % mod;
	return r;
}
// Large inputs: a(x) b(x) = c(x) at random points mod the prime `mod`,
// plus some coefficients by direct summation.
void checkLarge(const vl& a, const vl& b) {
	vl c = convMod<mod>(a, b);
	assert(sz(c) == sz(a) + sz(b) - 1);
	for (ll x : c) assert(0 <= x && x < mod);
	rep(it,0,4) {
		ll x = rnd(1, mod - 1);
		assert(mulM(evalM(a, x), evalM(b, x)) == evalM(c, x));
	}
	rep(it,0,20) {
		int i = (int)rnd(0, sz(c) - 1);
		ll s = 0;
		rep(j,max(0, i - sz(b) + 1),min(i, sz(a) - 1) + 1) s = (s + a[j] * b[i - j]) % mod;
		assert(s == c[i]);
	}
}

int main() {
	vl a, b;
	rep(it,0,6000) {
		a.resize(ra() % 100);
		b.resize(ra() % 100);
		for(auto &x: a) x = ra() % mod;
		for(auto &x: b) x = ra() % mod;
		auto v1 = simpleConv<mod>(a, b);
		auto v2 = convMod<mod>(a, b);
		assert(v1 == v2);
	}

	// Other moduli: tiny, non-prime, perfect squares, largest int.
	testSmall<1>(200, 20);
	testSmall<2>(2000, 40);
	testSmall<3>(2000, 40);
	testSmall<4>(2000, 40);
	testSmall<7>(2000, 40);
	testSmall<1000>(2000, 40);
	testSmall<1000000>(2000, 40);
	testSmall<998244353>(3000, 60);
	testSmall<1000000007>(3000, 60);
	testSmall<2147395600>(3000, 60); // 46340^2
	testSmall<2147483647>(3000, 60);

	// Worst case at the documented bound N log2(N) mod < 8.6e14:
	// N = |A| + |B| = 2^15 for mod = 1e9+7.
	{
		int n = 1 << 14;
		assert(2.0 * n * 15 * (double)mod < 8.6e14);
		a.assign(n, mod - 1); b.assign(n, mod - 1);
		checkLarge(a, b);
		for (auto& x : a) x = rnd(0, mod - 1);
		for (auto& x : b) x = rnd(0, mod - 1);
		checkLarge(a, b);
		for (auto& x : a) x = rnd(0, 1) * (mod - 1);
		for (auto& x : b) x = rnd(0, 1) * (mod - 1);
		checkLarge(a, b);
	}
	// Beyond the proven bound but within the "in practice" one (1e16):
	// N = 2^18.
	{
		int n = 1 << 17;
		a.assign(n, mod - 1); b.assign(n, mod - 1);
		checkLarge(a, b);
		for (auto& x : a) x = rnd(0, mod - 1);
		for (auto& x : b) x = rnd(0, mod - 1);
		checkLarge(a, b);
	}
	// Unbalanced, non power of two sizes.
	rep(it,0,30) {
		a.resize(rnd(1, 20000)); b.resize(rnd(1, 20000));
		for (auto& x : a) x = rnd(0, mod - 1);
		for (auto& x : b) x = rnd(0, mod - 1);
		checkLarge(a, b);
	}
	cout<<"Tests passed!"<<endl;
}
