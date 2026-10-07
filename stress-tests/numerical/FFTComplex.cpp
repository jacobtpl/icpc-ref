#include "../utilities/template.h"

// FFTComplex.h uses `C` and `fft` from FastFourierTransform.h, which it
// has to pull in itself.
#include "../../content/numerical/FFTComplex.h"

mt19937_64 rng(4242);
ll rnd(ll lo, ll hi) { return lo + (ll)(rng() % (unsigned long long)(hi - lo + 1)); }

vector<C> naive(const vector<C>& a, const vector<C>& b) {
	if (a.empty() || b.empty()) return {};
	vector<C> c(sz(a) + sz(b) - 1);
	rep(i,0,sz(a)) rep(j,0,sz(b)) c[i + j] += a[i] * b[j];
	return c;
}

const ll P = (1LL << 61) - 1;
ll mulP(ll a, ll b) { return (ll)((__int128)a * b % P); }
struct G { ll re, im; }; // Gaussian integers mod P
G mulG(G a, G b) {
	return {(mulP(a.re, b.re) - mulP(a.im, b.im) + P) % P,
		(mulP(a.re, b.im) + mulP(a.im, b.re)) % P};
}
G evalG(const vector<C>& v, ll x) {
	G r{0, 0};
	for (int i = sz(v) - 1; i >= 0; i--) {
		r.re = (mulP(r.re, x) + llround(v[i].real()) % P + P) % P;
		r.im = (mulP(r.im, x) + llround(v[i].imag()) % P + P) % P;
	}
	return r;
}

int main() {
	assert(conv_complex({}, {}).empty());
	assert(conv_complex({C(1, 2)}, {}).empty());
	assert(conv_complex({}, {C(1, 2)}).empty());
	{
		auto r = conv_complex({C(1, 2)}, {C(3, -1)});
		assert(sz(r) == 1 && abs(r[0] - C(5, 5)) < 1e-12);
	}
	// All small size combinations, complex values.
	rep(it,0,30000) {
		vector<C> a(rnd(1, 40)), b(rnd(1, 40));
		for (C& x : a) x = C((double)rnd(-5000000, 5000000) / 1e6, (double)rnd(-5000000, 5000000) / 1e6);
		for (C& x : b) x = C((double)rnd(-5000000, 5000000) / 1e6, (double)rnd(-5000000, 5000000) / 1e6);
		if (rnd(0, 4) == 0) for (C& x : a) x = x.real(); // purely real
		if (rnd(0, 4) == 0) for (C& x : b) x = C(0, x.imag()); // purely imaginary
		auto c = conv_complex(a, b), d = naive(a, b);
		assert(sz(c) == sz(d));
		rep(i,0,sz(c)) assert(abs(c[i] - d[i]) < 1e-8);
	}
	// Medium sizes, including results whose size is a power of two.
	for (int n : {127, 128, 129, 255, 256, 257, 1000, 2048}) {
		for (int m : {1, 2, n, n + 1}) {
			vector<C> a(n), b(m);
			for (C& x : a) x = C((double)rnd(-1000, 1000), (double)rnd(-1000, 1000));
			for (C& x : b) x = C((double)rnd(-1000, 1000), (double)rnd(-1000, 1000));
			auto c = conv_complex(a, b), d = naive(a, b);
			assert(sz(c) == sz(d));
			rep(i,0,sz(c)) assert(abs(c[i] - d[i]) < 1e-3);
		}
	}
	// Large: Gaussian integers, checked exactly after rounding by
	// evaluating a(x) b(x) = c(x) at random points mod 2^61-1.
	{
		int n = 1 << 18;
		vector<C> a(n), b(n);
		rep(type,0,2) {
			for (C& x : a) x = type ? C(1000, -1000) : C((double)rnd(-1000, 1000), (double)rnd(-1000, 1000));
			for (C& x : b) x = type ? C(-1000, 1000) : C((double)rnd(-1000, 1000), (double)rnd(-1000, 1000));
			auto c = conv_complex(a, b);
			assert(sz(c) == 2 * n - 1);
			for (C& x : c) {
				assert(abs(x - C((double)llround(x.real()), (double)llround(x.imag()))) < 0.2);
			}
			rep(it,0,3) {
				ll x = rnd(1, P - 1);
				G l = mulG(evalG(a, x), evalG(b, x)), r = evalG(c, x);
				assert(l.re == r.re && l.im == r.im);
			}
		}
	}
	cout<<"Tests passed!"<<endl;
}
