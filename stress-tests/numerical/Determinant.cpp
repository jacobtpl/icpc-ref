#include "../utilities/template.h"

#include "../../content/numerical/Determinant.h"

mt19937_64 rng(2024);
int rnd(int lo, int hi) { return lo + (int)(rng() % (unsigned)(hi - lo + 1)); }

// Exact oracle: fraction-free (Bareiss) elimination on integers.
ll exactDet(vector<vector<ll>> a) {
	int n = sz(a);
	ll sign = 1, prev = 1;
	rep(i,0,n) {
		int p = -1;
		rep(j,i,n) if (a[j][i]) { p = j; break; }
		if (p == -1) return 0;
		if (p != i) swap(a[p], a[i]), sign = -sign;
		rep(j,i+1,n) rep(k,i+1,n)
			a[j][k] = (ll)(((__int128)a[j][k] * a[i][i] - (__int128)a[j][i] * a[i][k]) / prev);
		prev = a[i][i];
	}
	return n ? sign * a[n-1][n-1] : 1;
}

// Independent floating oracle: Laplace expansion over subsets, O(2^n n).
long double laplaceDet(const vector<vector<double>>& a) {
	int n = sz(a);
	vector<long double> dp(1 << n);
	dp[0] = 1;
	rep(mask,0,1<<n) {
		int r = __builtin_popcount(mask);
		if (r == n) break;
		int sgn = 1;
		for (int c = n - 1; c >= 0; c--) {
			if (mask >> c & 1) sgn = -sgn;
			else dp[mask | 1 << c] += sgn * dp[mask] * a[r][c];
		}
	}
	return dp[(1 << n) - 1];
}

vector<vector<double>> toD(const vector<vector<ll>>& a) {
	vector<vector<double>> d(sz(a));
	rep(i,0,sz(a)) d[i].assign(all(a[i]));
	return d;
}

// A singular matrix does not give exactly 0 in floating point, so errors
// are measured relative to the Hadamard-style bound prod_i |row_i|_1.
double scale(const vector<vector<double>>& a) {
	double s = 1;
	for (auto& r : a) {
		double t = 0;
		for (double x : r) t += fabs(x);
		s *= t;
	}
	return max(s, 1.0);
}

int main() {
	{ // edge cases
		vector<vector<double>> e;
		assert(det(e) == 1);
		vector<vector<double>> one{{-3.5}};
		assert(det(one) == -3.5);
		vector<vector<double>> z(4, vector<double>(4));
		assert(det(z) == 0);
		vector<vector<double>> sw{{0, 1}, {1, 0}};
		assert(det(sw) == -1);
	}
	// Small integer matrices against the exact determinant.
	rep(it,0,300000) {
		int n = rnd(0, 7), lim = rnd(1, 10), type = rnd(0, 3);
		vector<vector<ll>> a(n, vector<ll>(n));
		for (auto& r : a) for (ll& x : r) x = rnd(-lim, lim);
		if (type == 1) for (auto& r : a) for (ll& x : r) if (rnd(0, 2)) x = 0; // sparse
		if (type == 2 && n >= 2) { // row = combination of two others
			int i = rnd(0, n-1), j = rnd(0, n-1), k = rnd(0, n-1);
			int c1 = rnd(-2, 2), c2 = rnd(-2, 2);
			if (i != j && i != k) rep(c,0,n) a[i][c] = c1 * a[j][c] + c2 * a[k][c];
		}
		if (type == 3 && n >= 2) { // duplicate column
			int i = rnd(0, n-1), j = rnd(0, n-1);
			rep(r,0,n) a[r][i] = a[r][j];
		}
		ll want = exactDet(a);
		auto d = toD(a), d2 = d;
		double got = det(d2);
		assert(fabs(got - (double)want) <= 1e-11 * scale(d));
	}
	// Integers near the limit of exactness in doubles.
	rep(it,0,20000) {
		int n = rnd(1, 4);
		vector<vector<ll>> a(n, vector<ll>(n));
		for (auto& r : a) for (ll& x : r) x = rnd(-30000, 30000);
		ll want = exactDet(a);
		auto d = toD(a), d2 = d;
		double got = det(d2);
		assert(fabs(got - (double)want) <= 1e-11 * scale(d));
	}
	// Real matrices against Laplace expansion in long double.
	rep(it,0,3000) {
		int n = rnd(1, 10);
		vector<vector<double>> a(n, vector<double>(n));
		for (auto& r : a) for (double& x : r) x = (double)(rng() % 2000001) / 1e6 - 1;
		long double want = laplaceDet(a);
		auto b = a;
		double got = det(b);
		assert(fabsl(got - want) <= 1e-11 * scale(a));
	}
	// Scaling: det(c * A) = c^n det(A), with huge and tiny magnitudes.
	rep(it,0,2000) {
		int n = rnd(1, 6);
		vector<vector<ll>> a(n, vector<ll>(n));
		for (auto& r : a) for (ll& x : r) x = rnd(-5, 5);
		double c = rnd(0, 1) ? 1e30 : 1e-30;
		auto d = toD(a);
		for (auto& r : d) for (double& x : r) x *= c;
		double want = (double)exactDet(a) * pow(c, n), got = det(d);
		assert(fabs(got - want) <= 1e-11 * scale(toD(a)) * pow(c, n));
	}
	cout<<"Tests passed!"<<endl;
}
