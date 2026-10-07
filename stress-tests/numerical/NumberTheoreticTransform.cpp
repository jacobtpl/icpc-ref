#include "../utilities/template.h"

typedef vector<ll> vl;
namespace ignore {
#include "../../content/number-theory/ModPow.h"
}
ll modpow(ll a, ll e);
#include "../../content/numerical/NumberTheoreticTransform.h"
ll modpow(ll a, ll e) {
    if (e == 0)
        return 1;
    ll x = modpow(a * a % mod, e >> 1);
    return e & 1 ? x * a % mod : x;
}

vl simpleConv(vl a, vl b) {
	int s = sz(a) + sz(b) - 1;
	if (a.empty() || b.empty()) return {};
	vl c(s);
	rep(i,0,sz(a)) rep(j,0,sz(b))
		c[i+j] = (c[i+j] + (ll)a[i] * b[j]) % mod;
	for(auto &x: c) if (x < 0) x += mod;
	return c;
}

int ra() {
    static unsigned X;
    X *= 123671231;
    X += 1238713;
    X ^= 1237618;
    return (X >> 1);
}

int main() {
	ll res = 0, res2 = 0;
	int ind = 0, ind2 = 0;
	vl a, b;
	rep(it,0,6000) {
		a.resize(ra() % 10);
		b.resize(ra() % 10);
		for(auto &x: a) x = (ra() % 100 - 50+mod)%mod;
		for(auto &x: b) x = (ra() % 100 - 50+mod)%mod;
		for(auto &x: simpleConv(a, b)) res += (ll)x * ind++ % mod;
		for(auto &x: conv(a, b)) res2 += (ll)x * ind2++ % mod;
		a.resize(16);
        vl a2 = a;
        ntt(a2);
        rep(k, 0, sz(a2)) {
            ll sum = 0;
            rep(x, 0, sz(a2)) { sum = (sum + a[x] * modpow(root, k * x * (mod - 1) / sz(a))) % mod; }
            assert(sum == a2[k]);
        }
	}
	assert(res==res2);
	// larger sizes, full-range and extreme values, power-of-two boundaries
	mt19937_64 rng(5);
	rep(it,0,300) {
		int n = (int)(rng() % 600) + 1, m = (int)(rng() % 600) + 1, mode = (int)(rng() % 3);
		if (it < 40) n = 1 << (it % 10), m = (1 << (it % 10)) + it / 10 % 3 - 1;
		if (m < 1) m = 1;
		a.resize(n), b.resize(m);
		for(auto &x: a) x = mode == 0 ? mod - 1 : mode == 1 ? (ll)(rng() % mod) : (rng() % 4 ? 0 : mod - 1);
		for(auto &x: b) x = mode == 0 ? mod - 1 : (ll)(rng() % mod);
		vl c = conv(a, b);
		for (auto x : c) assert(0 <= x && x < mod);
		assert(c == simpleConv(a, b));
	}
	assert(conv({}, {1, 2}).empty() && conv({1, 2}, {}).empty());
	// transform sizes in non-monotone order (the root table is static), forward then inverse
	for (int n : {1, 2, 1024, 4, 1 << 16, 8, 1 << 12, 1}) {
		vl v(n);
		for (auto &x : v) x = (ll)(rng() % mod);
		vl w = v;
		ntt(w);
		if (n <= 1024) rep(k,0,n) {
			ll sum = 0, g = modpow(root, (mod - 1) / n * k), pw = 1;
			rep(x,0,n) sum = (sum + v[x] * pw) % mod, pw = pw * g % mod;
			assert(sum == w[k]);
		}
		reverse(w.begin() + 1, w.end());
		ntt(w);
		ll in = modpow(n, mod - 2);
		rep(i,0,n) assert(w[i] * in % mod == v[i]);
	}
	{ // 2^19 x 2^19, checked by evaluating at random points
		int n = 1 << 19;
		a.resize(n), b.resize(n - 3);
		for(auto &x: a) x = (ll)(rng() % mod);
		for(auto &x: b) x = rng() % 2 ? mod - 1 : (ll)(rng() % mod);
		vl c = conv(a, b);
		assert(sz(c) == sz(a) + sz(b) - 1);
		rep(it,0,3) {
			ll x = (ll)(rng() % mod);
			auto ev = [&](const vl& p) { ll r = 0; for (int i = sz(p); i--;) r = (r * x + p[i]) % mod; return r; };
			assert(ev(a) * ev(b) % mod == ev(c));
		}
	}
	cout<<"Tests passed!"<<endl;
}
