#include "../utilities/template.h"

#include "../../content/numerical/PolyRoots.h"

mt19937_64 rng(11);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }
double rndd(double lo, double hi) { return uniform_real_distribution<double>(lo, hi)(rng); }

Poly mulLin(const Poly& p, double r) { // p * (x - r)
	Poly q; q.a.assign(sz(p.a) + 1, 0);
	rep(i,0,sz(p.a)) q.a[i + 1] += p.a[i], q.a[i] -= r * p.a[i];
	return q;
}
Poly mulQuad(const Poly& p, double b, double c) { // p * (x^2 + bx + c)
	Poly q; q.a.assign(sz(p.a) + 2, 0);
	rep(i,0,sz(p.a)) q.a[i + 2] += p.a[i], q.a[i + 1] += b * p.a[i], q.a[i] += c * p.a[i];
	return q;
}

int main() {
	// Usage example from the header
	{
		auto r = polyRoots({{2, -3, 1}}, -1e9, 1e9);
		assert(sz(r) == 2 && fabs(r[0] - 1) < 1e-9 && fabs(r[1] - 2) < 1e-9);
	}
	// distinct well-separated real roots (plus optional complex pairs), all inside the range
	rep(it,0,40000) {
		int k = (int)rnd(1, 6), cp = (int)rnd(0, 2);
		vector<double> roots;
		bool integral = rnd(0, 1);
		while (sz(roots) < k) {
			double r = integral ? (double)rnd(-10, 10) : rndd(-10, 10);
			bool ok = 1;
			for (double o : roots) if (fabs(o - r) < 0.25) ok = 0;
			if (ok) roots.push_back(r);
		}
		sort(all(roots));
		Poly p{{(double)rnd(1, 5) * (rnd(0, 1) ? 1 : -1)}};
		for (double r : roots) p = mulLin(p, r);
		rep(i,0,cp) { double b = rndd(-3, 3); p = mulQuad(p, b, b * b / 4 + rndd(0.5, 3)); }
		double lo = rnd(0, 1) ? -1e9 : -11, hi = rnd(0, 1) ? 1e9 : 11;
		auto got = polyRoots(p, lo, hi);
		assert(is_sorted(all(got)));
		if (sz(got) != k) { cerr << "expected " << k << " roots, got " << sz(got) << endl; assert(0); }
		rep(i,0,k) assert(fabs(got[i] - roots[i]) < 1e-6);
	}
	// no real roots
	rep(it,0,2000) {
		Poly p{{1}};
		rep(i,0,(int)rnd(1, 3)) { double b = rndd(-3, 3); p = mulQuad(p, b, b * b / 4 + rndd(0.5, 3)); }
		assert(polyRoots(p, -1e9, 1e9).empty());
	}
	// a root inside [xmin, xmax] is always found, even when the other roots are far outside
	rep(it,0,5000) {
		double a = (double)rnd(-20, -5), b = (double)rnd(-2, 2), c = (double)rnd(5, 20);
		Poly p = mulLin(mulLin(mulLin(Poly{{1}}, a), b), c);
		auto got = polyRoots(p, -4, 4);
		bool found = 0;
		for (double r : got) found |= fabs(r - b) < 1e-7;
		assert(found);
	}
#ifdef POLYROOTS_STRICT_RANGE
	// Roots outside [xmin, xmax] are reported: here both roots (1 and 2) lie outside [5, 10].
	assert(polyRoots({{2, -3, 1}}, 5, 10).empty());
	assert(polyRoots({{-1, 1}}, 5, 10).empty());
#endif
#ifdef POLYROOTS_MULTIPLE
	// Roots of even multiplicity produce no sign change and are missed.
	assert(sz(polyRoots({{1, -2, 1}}, -10, 10)) == 1);
#endif
	cout << "Tests passed!" << endl;
}
