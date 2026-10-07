#include "../utilities/template.h"

#include "../../content/numerical/PolyInterpolate.h"

mt19937_64 rng(13);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }

int main() {
	// small integer polynomials through distinct integer points: compare coefficients
	rep(it,0,100000) {
		int n = (int)rnd(1, 7);
		vector<ll> c(n);
		for (auto& v : c) v = rnd(-20, 20);
		vi xs;
		while (sz(xs) < n) { int v = (int)rnd(-6, 6); if (!count(all(xs), v)) xs.push_back(v); }
		vd x(n), y(n);
		rep(i,0,n) {
			ll val = 0;
			for (int j = n; j--;) val = val * xs[i] + c[j];
			x[i] = xs[i], y[i] = (double)val;
		}
		vd res = interpolate(x, y, n);
		assert(sz(res) == n);
		rep(i,0,n) assert(fabs(res[i] - (double)c[i]) < 1e-6);
	}
	// Chebyshev nodes as recommended in the header: result must reproduce the points
	// (coefficient error grows about 4x per extra point: ~1e-8 at n=14, ~1e-4 at n=20)
	for (int n = 2; n <= 14; n++) rep(it,0,300) {
		vd x(n), y(n), c(n);
		for (auto& v : c) v = (double)rnd(-1000, 1000) / 100;
		rep(k,0,n) {
			x[k] = cos(k * acos(-1.0) / (n - 1));
			double val = 0;
			for (int j = n; j--;) val = val * x[k] + c[j];
			y[k] = val;
		}
		vd res = interpolate(x, y, n);
		rep(i,0,n) assert(fabs(res[i] - c[i]) < 1e-6);
	}
	{ // n = 1
		vd res = interpolate({5}, {-3}, 1);
		assert(sz(res) == 1 && res[0] == -3);
	}
	cout << "Tests passed!" << endl;
}
