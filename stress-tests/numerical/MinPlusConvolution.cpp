#include "../utilities/template.h"

#include "../../content/numerical/MinPlusConvolution.h"

vl naive(const vl& a, const vl& b) {
	if (a.empty() || b.empty()) return {};
	vl c(sz(a) + sz(b) - 1, LLONG_MAX);
	rep(i,0,sz(a)) rep(j,0,sz(b)) c[i+j] = min(c[i+j], a[i] + b[j]);
	return c;
}

ll rnd(ll lo, ll hi) { return lo + (ll)(rand() % (hi - lo + 1)); }

vl convex(int n, ll vlo, ll vhi, ll slo, ll shi) {
	if (!n) return {};
	vl d(n - 1);
	for (ll& x : d) x = rnd(slo, shi);
	sort(all(d));
	vl a(n);
	a[0] = rnd(vlo, vhi);
	rep(i,1,n) a[i] = a[i-1] + d[i-1];
	return a;
}

vl arbitrary(int n, ll vlo, ll vhi) {
	vl b(n);
	for (ll& x : b) x = rnd(vlo, vhi);
	return b;
}

void check(const vl& a, const vl& b, bool bConvex) {
	vl ex = naive(a, b);
	assert(minPlus(a, b) == ex);
	if (bConvex) {
		assert(minPlusCC(a, b) == ex);
		assert(minPlusCC(b, a) == ex);
		assert(minPlus(b, a) == ex);
	}
}

int main() {
	// Edge cases: empty inputs, single elements
	assert(minPlusCC({}, {}).empty());
	assert(minPlus({}, {}).empty());
	assert(minPlusCC({1, 2}, {}).empty());
	assert(minPlus({}, {5}).empty());
	assert(minPlusCC({3}, {4}) == vl{7});
	assert(minPlus({3}, {4}) == vl{7});
	assert(minPlus({3}, {4, -1, 2}) == (vl{7, 2, 5}));

	// Small random cases, heavy ties: slopes in {-1,0,1}, small values
	rep(it,0,30000) {
		int n = rand() % 8, m = rand() % 8;
		vl a = convex(n, -3, 3, -1, 1);
		check(a, convex(m, -3, 3, -1, 1), true);
		check(a, arbitrary(m, -3, 3), false);
	}
	// Medium random cases
	rep(it,0,3000) {
		int n = rand() % 60, m = rand() % 60;
		vl a = convex(n, -100, 100, -20, 20);
		check(a, convex(m, -100, 100, -20, 20), true);
		check(a, arbitrary(m, -1000, 1000), false);
	}
	// Large values: every sum is close to the ll limits
	const ll B = (ll)4.5e18, S = (ll)1e17;
	rep(it,0,3000) {
		int n = rand() % 30, m = rand() % 30;
		vl a = convex(n, -B + 40 * S, B - 40 * S, -S, S);
		if (n > 1) rep(i,0,n) assert(-B <= a[i] && a[i] <= B);
		check(a, convex(m, -B + 40 * S, B - 40 * S, -S, S), true);
		check(a, arbitrary(m, -B, B), false);
		check(a, arbitrary(m, B - 1000, B), false);
		check(a, arbitrary(m, -B, -B + 1000), false);
	}
	// Degenerate: all-zero / linear / all-equal b, very unbalanced sizes
	rep(it,0,300) {
		int n = rand() % 200 + 1, m = rand() % 3 + 1;
		if (rand() % 2) swap(n, m);
		vl z(n, 0), lin(n), cst(m, rnd(-5, 5));
		rep(i,0,n) lin[i] = 3 * i - 7;
		check(z, cst, true);
		check(lin, cst, true);
		check(z, arbitrary(m, -5, 5), false);
		check(lin, arbitrary(m, -5, 5), false);
		check(convex(n, -50, 50, -5, 5), arbitrary(m, -50, 50), false);
	}
	// A few larger cases against brute force
	rep(it,0,20) {
		int n = rand() % 1500 + 1, m = rand() % 1500 + 1;
		vl a = convex(n, -1000000, 1000000, -1000, 1000);
		check(a, convex(m, -1000000, 1000000, -1000, 1000), true);
		check(a, arbitrary(m, -1000000, 1000000), false);
	}
	cout<<"Tests passed!"<<endl;
}
