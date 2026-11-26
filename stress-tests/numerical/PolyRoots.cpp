#include "../utilities/template.h"
#include "../../content/numerical/PolyRoots.h"

void check(vector<double> got, vector<double> want) {
	sort(all(want));
	assert(sz(got) == sz(want));
	rep(i,0,sz(got)) assert(abs(got[i] - want[i]) < 1e-6);
}

int main() {
	mt19937 rng(777);
	rep(it,0,200000) {
		// Distinct roots at multiples of 0.5 in [-10, 10], plus optional
		// irreducible quadratic factors and leading zero coefficients.
		int k = (int)(rng() % 7);
		set<int> rs;
		while (sz(rs) < k) rs.insert((int)(rng() % 41) - 20);
		vector<double> a = {(rng() % 2 ? 1.0 : -1.0) * (1 + (int)(rng() % 5))};
		auto mul = [&](vector<double> f) {
			vector<double> c(sz(a) + sz(f) - 1);
			rep(i,0,sz(a)) rep(j,0,sz(f)) c[i+j] += a[i] * f[j];
			a = c;
		};
		vector<double> roots;
		for (int r : rs) roots.push_back(r / 2.0), mul({-r / 2.0, 1});
		rep(q,0,(int)(rng() % 2)) mul({1.0 + (int)(rng() % 3), 1.0 * (int)(rng() % 2), 1});
		rep(z,0,(int)(rng() % 3)) a.push_back(0);
		check(polyRoots({a}, -1e3, 1e3), roots);
		// A range with endpoints strictly between candidate roots.
		double lo = (int)(rng() % 41) / 2.0 - 10.25, hi = lo + (int)(rng() % 20);
		vector<double> in;
		for (double r : roots) if (lo < r && r < hi) in.push_back(r);
		check(polyRoots({a}, lo, hi), in);
	}
	// Roots exactly on the endpoints are included.
	check(polyRoots({{2,-3,1}}, 1, 2), {1, 2});
	check(polyRoots({{2,-3,1}}, 1, 1.5), {1});
	check(polyRoots({{2,-3,1}}, 1.5, 2), {2});
	check(polyRoots({{-6,11,-6,1}}, 1, 3), {1, 2, 3});
	check(polyRoots({{-6,11,-6,1}}, 2, 2), {2});
	check(polyRoots({{-5,1}}, 0, 1), {});
	check(polyRoots({{-5,1}}, 0, 5), {5});
	// Degenerate polynomials.
	check(polyRoots({{}}, -1, 1), {});
	check(polyRoots({{0}}, -1, 1), {});
	check(polyRoots({{3}}, -1, 1), {});
	check(polyRoots({{0,0,0}}, -1, 1), {});
	check(polyRoots({{1,0}}, -1, 1), {});
	cout << "Tests passed!" << endl;
}
