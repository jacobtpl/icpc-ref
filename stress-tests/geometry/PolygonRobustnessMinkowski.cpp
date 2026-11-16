#include "../utilities/template.h"

#include "../../content/geometry/Point.h"
#include "../../content/geometry/MinkowskiSum.h"

// minkowski_sum against the strict hull of all pairwise sums, on inputs
// with collinear and repeated boundary vertices, points and segments.
typedef Point<ll> P;
typedef __int128 lll;
mt19937_64 rng(20261013);
ll rnd(ll a, ll b) { return uniform_int_distribution<ll>(a, b)(rng); }

lll X(P a, P b, P c) {
	return (lll)(b.x - a.x) * (c.y - a.y) - (lll)(b.y - a.y) * (c.x - a.x);
}
// strict hull, CCW from the lexicographically smallest point
vector<P> hull(vector<P> v) {
	sort(all(v)); v.erase(unique(all(v)), v.end());
	if (sz(v) <= 2) return v;
	vector<P> h(sz(v)+1);
	int s = 0, t = 0;
	for (int it = 2; it--; s = --t, reverse(all(v)))
		for (P p : v) {
			while (t >= s + 2 && X(h[t-2], h[t-1], p) <= 0) t--;
			h[t++] = p;
		}
	h.resize(t);
	return h;
}
vector<P> brute(const vector<P>& a, const vector<P>& b) {
	vector<P> s;
	for (P p : a) for (P q : b) s.push_back(p + q);
	return hull(s);
}
// convex polygon (possibly a point or segment) with extra boundary points
vector<P> gen(int n, ll R, ll mul) {
	vector<P> v;
	int type = (int)rnd(0, 5);
	if (type == 0) v = {P(rnd(-R, R), rnd(-R, R))};
	else if (type == 1) { // collinear points
		P o(rnd(-R/2, R/2), rnd(-R/2, R/2)), d(rnd(-3, 3), rnd(-3, 3));
		ll k = max<ll>(R / 8, 1);
		rep(i,0,n) v.push_back(o + d * rnd(-k, k));
	} else rep(i,0,n) v.push_back(P(rnd(-R, R), rnd(-R, R)));
	vector<P> h = hull(v), res;
	int m = sz(h);
	rep(i,0,m) {
		P p = h[i], q = h[(i+1)%m];
		res.push_back(p);
		if (rnd(0, 2) == 0) res.push_back(p); // repeated vertex
		if (m < 2) continue;
		// lattice points strictly inside the edge, in order
		ll g = __gcd(abs(q.x - p.x), abs(q.y - p.y));
		P d((q.x - p.x) / g, (q.y - p.y) / g);
		vector<ll> ks;
		rep(j,0,(int)rnd(0, 2)) if (g > 1) ks.push_back(rnd(1, g - 1));
		sort(all(ks));
		for (ll k : ks) res.push_back(p + d * k);
	}
	rotate(res.begin(), res.begin() + rnd(0, sz(res) - 1), res.end());
	for (P& p : res) p = p * mul;
	return res;
}
void check(const vector<P>& a, const vector<P>& b) {
	vector<P> want = brute(a, b), got = minkowski_sum(a, b);
	if (got != want) {
		cerr << "minkowski mismatch\na:";
		for (P p : a) cerr << ' ' << p;
		cerr << "\nb:";
		for (P p : b) cerr << ' ' << p;
		cerr << "\nwant:";
		for (P p : want) cerr << ' ' << p;
		cerr << "\ngot:";
		for (P p : got) cerr << ' ' << p;
		cerr << endl;
		abort();
	}
	assert(minkowski_sum(b, a) == want);
}

int main() {
	// --- review reproducers ---
	// segment + square with a collinear boundary vertex (was out of bounds)
	check({P(0,0), P(2,0)}, {P(0,0), P(1,0), P(2,0), P(2,2), P(0,2)});
	check({P(0,0), P(1,0), P(2,0), P(2,2), P(0,2)},
		{P(0,0), P(1,0), P(2,0), P(2,2), P(0,2)});
	// all edges parallel in both inputs
	check({P(0,0), P(1,1), P(2,2)}, {P(0,0), P(3,3)});
	check({P(0,0), P(0,5)}, {P(0,0), P(0,1), P(0,7)});
	check({P(0,0), P(0,5)}, {P(0,0), P(4,0)});
	check({P(3,4)}, {P(-1,2)});
	check({P(3,4)}, {P(0,0), P(1,0), P(1,1), P(0,1)});
	check({P(3,4), P(3,4), P(3,4)}, {P(1,1), P(1,1)});
	// empty input
	assert(minkowski_sum(vector<P>{}, vector<P>{P(1,1)}).empty());
	assert(minkowski_sum(vector<P>{P(1,1)}, vector<P>{}).empty());
	assert(minkowski_sum(vector<P>{}, vector<P>{}).empty());
	// any starting vertex
	{
		vector<P> sq = {P(0,0), P(1,0), P(1,1), P(0,1)};
		vector<P> tr = {P(0,0), P(2,1), P(-1,3)};
		rep(i,0,4) rep(j,0,3) {
			check(sq, tr);
			rotate(tr.begin(), tr.begin() + 1, tr.end());
			if (j == 2) rotate(sq.begin(), sq.begin() + 1, sq.end());
		}
	}
	// doubles
	{
		typedef Point<double> D;
		vector<D> a = {D(0,0), D(0.5,0), D(1,0), D(1,1), D(0,1)}, b = {D(1e9,1e9), D(1e9+1,1e9)};
		vector<D> want = {D(1e9,1e9), D(1e9+2,1e9), D(1e9+2,1e9+1), D(1e9,1e9+1)};
		assert(minkowski_sum(a, b) == want);
	}
	// --- randomized ---
	ll cnt[4] = {};
	rep(it,0,300000) {
		ll R = it % 3 == 0 ? 2 : it % 3 == 1 ? 6 : 40;
		vector<P> a = gen((int)rnd(1, 8), R, 1), b = gen((int)rnd(1, 8), R, 1);
		check(a, b);
		cnt[min(sz(brute(a, b)), 3)]++;
	}
	rep(i,1,4) assert(cnt[i] > 1000);
	// coordinates up to 2^30 - 1 (documented exact range)
	rep(it,0,50000) {
		ll mul = ((1LL << 30) - 1) / 40;
		vector<P> a = gen((int)rnd(1, 8), 40, mul), b = gen((int)rnd(1, 8), 40, mul);
		check(a, b);
	}
	{
		ll M = (1LL << 30) - 1;
		vector<P> a = {P(-M,-M), P(0,-M), P(M,-M), P(M,M), P(-M,M)};
		check(a, a);
		check(a, {P(-M,M), P(M,-M)});
		check({P(-M,-M), P(M,M)}, {P(-M,M), P(M,-M)});
	}
	// --- tiny repeat and large input ---
	{
		vector<P> sq = {P(0,0), P(1,0), P(1,1), P(0,1)}, sg = {P(0,0), P(1,1)};
		ll s = 0;
		rep(it,0,300000) s += sz(minkowski_sum(sq, sg));
		assert(s == 6 * 300000LL);
	}
	{
		// 2*10^5 lattice points on a parabola-like convex chain
		auto mk = [](int n, ll a, ll b) {
			vector<P> up, dn;
			rep(i,0,n) {
				ll x = i - n / 2;
				dn.push_back(P(a * x, b * x * x));
				up.push_back(P(a * x, 2 * b * (ll)n * n - b * x * x));
			}
			reverse(all(up));
			dn.insert(dn.end(), all(up));
			return dn;
		};
		vector<P> a = mk(100000, 3, 1), b = mk(100000, 2, 1);
		auto r = minkowski_sum(a, b);
		assert(sz(r) > 300000 && sz(r) <= sz(a) + sz(b));
		rep(i,0,sz(r)) assert(X(r[i], r[(i+1)%sz(r)], r[(i+2)%sz(r)]) > 0);
		assert(r[0] == *min_element(all(a)) + *min_element(all(b)));
	}
	cout << "Tests passed!" << endl;
}
