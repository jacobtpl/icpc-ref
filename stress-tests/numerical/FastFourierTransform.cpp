#include "../utilities/template.h"
#include <valarray>

#include "../../content/numerical/FastFourierTransform.h"

const double eps = 1e-8;
mt19937_64 rng(777);
ll rnd(ll lo, ll hi) { return lo + (ll)(rng() % (unsigned long long)(hi - lo + 1)); }

const ll P = (1LL << 61) - 1;
ll mulP(ll a, ll b) { return (ll)((__int128)a * b % P); }
// Value of the integer polynomial v at x, modulo the prime 2^61-1.
ll evalP(const vector<ll>& v, ll x) {
	ll r = 0;
	for (int i = sz(v) - 1; i >= 0; i--) r = (mulP(r, x) + v[i] % P + P) % P;
	return r;
}

// Checks round(conv(a, b)) exactly: a(x) b(x) = c(x) at random points
// mod P (Schwartz-Zippel), plus some coefficients by direct summation.
void checkExact(const vector<ll>& a, const vector<ll>& b) {
	vd c = conv(vd(all(a)), vd(all(b)));
	assert(sz(c) == sz(a) + sz(b) - 1);
	vector<ll> ci(sz(c));
	double worst = 0;
	rep(i,0,sz(c)) {
		ci[i] = llround(c[i]);
		worst = max(worst, fabs(c[i] - (double)ci[i]));
	}
	assert(worst < 0.3);
	rep(it,0,3) {
		ll x = rnd(1, P - 1);
		assert(mulP(evalP(a, x), evalP(b, x)) == evalP(ci, x));
	}
	rep(it,0,20) {
		int i = (int)rnd(0, sz(c) - 1);
		ll s = 0;
		rep(j,max(0, i - sz(b) + 1),min(i, sz(a) - 1) + 1) s += a[j] * b[i - j];
		assert(s == ci[i]);
	}
}

int main() {
	// fft against the definition, for all sizes in mixed order (the root
	// table is static and grows lazily).
	rep(it,0,3000) {
		int n = 1 << rnd(0, 8);
		vector<C> a(n);
		rep(i,0,n) a[i] = C((double)rnd(-5, 5), (double)rnd(-5, 5));
		auto aorig = a;
		fft(a);
		rep(k,0,n) {
			C sum{};
			rep(x,0,n) sum += aorig[x] * polar(1.0, 2 * M_PI * k * x / n);
			assert(norm(sum - a[k]) < 1e-12 * n * n * 50);
		}
	}

	// conv against the naive product: empty inputs, single elements,
	// all size combinations, real values.
	assert(conv({}, {}).empty() && conv({1}, {}).empty() && conv({}, {1}).empty());
	rep(it,0,20000) {
		vd A(rnd(1, 40)), B(rnd(1, 40));
		for(auto &x: A) x = (double)rnd(-5000000, 5000000) / 1e6;
		for(auto &x: B) x = (double)rnd(-5000000, 5000000) / 1e6;
		vd C = conv(A, B);
		assert(sz(C) == sz(A) + sz(B) - 1);
		rep(i,0,sz(A) + sz(B) - 1) {
			double sum = 0;
			rep(j,0,sz(A)) if (i - j >= 0 && i - j < sz(B)) {
				sum += A[j] * B[i - j];
			}
			assert(abs(sum - C[i]) < eps);
		}
	}

	// Rounding at the documented bound
	// (sum a_i^2 + sum b_i^2) log2(N) < 9e14 with N = |A| + |B| = 2^20:
	// that allows |a_i|, |b_i| <= 6500.
	{
		int n = 1 << 19, v = 6500;
		assert(2.0 * n * v * v * 20 < 9e14);
		vector<ll> a(n), b(n);
		for (ll& x : a) x = rnd(-v, v);
		for (ll& x : b) x = rnd(-v, v);
		checkExact(a, b); // random
		a.assign(n, v); b.assign(n, v);
		checkExact(a, b); // all maximal
		for (ll& x : a) x = rnd(0, 1) ? v : -v;
		for (ll& x : b) x = rnd(0, 1) ? v : -v;
		checkExact(a, b); // random signs
		rep(i,0,n) a[i] = i % 2 ? v : -v, b[i] = v;
		checkExact(a, b); // alternating
	}
	// Unbalanced and non power of two sizes, non-negative values.
	rep(it,0,30) {
		vector<ll> a(rnd(1, 30000)), b(rnd(1, 30000));
		for (ll& x : a) x = rnd(0, 30000);
		for (ll& x : b) x = rnd(0, 30000);
		checkExact(a, b);
	}
	cout<<"Tests passed!"<<endl;
}
