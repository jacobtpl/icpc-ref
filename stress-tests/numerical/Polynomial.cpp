#include "../utilities/template.h"

#include "../../content/numerical/Polynomial.h"

mt19937_64 rng(7);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }

int main() {
	// Integer coefficients/points small enough that every double operation is exact.
	rep(it,0,200000) {
		int n = (int)rnd(1, 9);
		vector<ll> c(n);
		for (auto& x : c) x = rnd(-9, 9);
		if (rnd(0, 3) == 0) c.back() = 0; // leading zeros are allowed
		Poly p; p.a.assign(all(c));
		ll x0 = rnd(-6, 6);
		ll val = 0, pw = 1;
		rep(i,0,n) val += c[i] * pw, pw *= x0;
		assert(p(double(x0)) == double(val));

		Poly d = p; d.diff();
		assert(sz(d.a) == n - 1);
		rep(i,1,n) assert(d.a[i - 1] == double(i * c[i]));

		// p(x) = (x - x0) q(x) + p(x0)
		Poly q = p; q.divroot(double(x0));
		assert(sz(q.a) == n - 1);
		vector<ll> back(n);
		rep(i,0,n-1) {
			ll qi = (ll)q.a[i]; assert(double(qi) == q.a[i]);
			back[i + 1] += qi, back[i] -= x0 * qi;
		}
		back[0] += val;
		assert(back == c);
	}
	// non-integer evaluation against long double Horner
	rep(it,0,100000) {
		int n = (int)rnd(0, 12);
		Poly p; p.a.resize(n);
		for (auto& x : p.a) x = (double)rnd(-1000000, 1000000) / 1000.0;
		double x = (double)rnd(-3000, 3000) / 1000.0;
		long double v = 0, mag = 0;
		for (int i = n; i--;) v = v * x + p.a[i], mag = mag * fabsl(x) + fabsl(p.a[i]);
		assert(fabsl(p(x) - v) <= 1e-13L * mag + 1e-300L);
	}
	assert(Poly{}(3.0) == 0);
	cout << "Tests passed!" << endl;
}
