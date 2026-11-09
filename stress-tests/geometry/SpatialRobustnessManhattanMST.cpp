// Regression tests for ManhattanMST.h: 64-bit weights (the old int version
// returned -294967296 for (-1e9,-1e9),(1e9,1e9)), extreme int/ll coordinates,
// int and ll inputs, degenerate inputs and many tiny testcases.
#include "../utilities/template.h"

#include "../../content/geometry/Point.h"
#include "../../content/geometry/ManhattanMST.h"
#include "../../content/data-structures/UnionFind.h"

typedef __int128 lll;
mt19937_64 rng(20261012);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }

template<class T> lll manh(Point<T> a, Point<T> b) {
	lll dx = (lll)a.x - b.x, dy = (lll)a.y - b.y;
	return (dx < 0 ? -dx : dx) + (dy < 0 ? -dy : dy);
}

// O(n^2) Prim with 128-bit distances
template<class T> lll bruteMST(const vector<Point<T>>& ps) {
	int n = sz(ps);
	if (n == 0) return 0;
	const lll BIG = (lll)1 << 100;
	vector<lll> d(n, BIG);
	vector<bool> in(n);
	d[0] = 0;
	lll tot = 0;
	rep(it,0,n) {
		int b = -1;
		rep(i,0,n) if (!in[i] && (b == -1 || d[i] < d[b])) b = i;
		in[b] = 1; tot += d[b];
		rep(i,0,n) if (!in[i]) d[i] = min(d[i], manh(ps[b], ps[i]));
	}
	return tot;
}

template<class T> void check(const vector<Point<T>>& ps) {
	int n = sz(ps);
	vector<array<ll, 3>> edges = manhattanMST(ps);
	assert(edges.size() <= 4 * ps.size());
	for (auto e : edges) {
		assert(0 <= e[1] && e[1] < n && 0 <= e[2] && e[2] < n);
		assert(e[1] != e[2]);
		assert(e[0] == manh(ps[e[1]], ps[e[2]]));
	}
	sort(all(edges));
	UF uf(n);
	lll cost = 0;
	int joined = 0;
	for (auto e : edges) if (uf.join((int)e[1], (int)e[2])) cost += e[0], joined++;
	assert(joined == max(n - 1, 0));
	assert(cost == bruteMST(ps));
}

template<class T> vector<Point<T>> gen(int n, ll lo, ll hi, int mode) {
	vector<Point<T>> ps;
	rep(i,0,n) {
		ll x = rnd(lo, hi), y = rnd(lo, hi);
		if (mode == 1) { // extremes only: many ties and maximal distances
			x = rnd(0, 1) ? lo : hi; y = rnd(0, 1) ? lo : hi;
			if (rnd(0, 3) == 0) x = rnd(lo, hi);
		}
		if (mode == 2) y = x;                     // diagonal
		if (mode == 3) y = lo + (hi - x);         // anti-diagonal
		if (mode == 4) y = lo;                    // horizontal line
		if (mode == 5 && i && rnd(0, 1)) {        // duplicates
			auto q = ps[rnd(0, i - 1)]; x = q.x; y = q.y;
		}
		ps.push_back(Point<T>((T)x, (T)y));
	}
	return ps;
}

// every candidate edge of a two-point input must carry the true distance
bool allWeights(const vector<array<ll, 3>>& e, ll w) {
	for (auto x : e) if (x[0] != w) return false;
	return !e.empty() && sz(e) <= 8;
}

int main() {
	const int IMIN = numeric_limits<int>::min(), IMAX = numeric_limits<int>::max();
	const ll LIM = 2'000'000'000'000'000'000LL; // documented ll bound

	{ // the review reproducer: true weight 4e9
		vector<Point<int>> a{Point<int>(-1000000000, -1000000000),
			Point<int>(1000000000, 1000000000)};
		auto e = manhattanMST(a);
		assert(allWeights(e, 4000000000LL));
		vector<Point<ll>> b{Point<ll>(-1000000000, -1000000000),
			Point<ll>(1000000000, 1000000000)};
		assert(manhattanMST(b) == e);
	}
	{ // full int range, including INT_MIN (which cannot be negated in int)
		vector<Point<int>> a{Point<int>(IMIN, IMIN), Point<int>(IMAX, IMAX),
			Point<int>(IMIN, IMAX), Point<int>(IMAX, IMIN), Point<int>(0, IMIN)};
		check(a);
		vector<Point<int>> b{Point<int>(IMIN, IMIN), Point<int>(IMAX, IMAX)};
		auto e = manhattanMST(b);
		assert(allWeights(e, 2 * (2LL * IMAX + 1)));
	}
	{ // documented ll bound: weight 8e18 still fits
		vector<Point<ll>> a{Point<ll>(-LIM, -LIM), Point<ll>(LIM, LIM)};
		auto e = manhattanMST(a);
		assert(allWeights(e, 4 * LIM));
		vector<Point<ll>> b{Point<ll>(-LIM, LIM), Point<ll>(LIM, -LIM),
			Point<ll>(LIM, LIM), Point<ll>(-LIM, -LIM), Point<ll>(0, 0)};
		check(b);
	}
	{ // empty, single, all-equal
		check(vector<Point<int>>{});
		check(vector<Point<ll>>{});
		assert(manhattanMST(vector<Point<int>>{Point<int>(IMIN, IMAX)}).empty());
		check(vector<Point<ll>>(5, Point<ll>(-LIM, LIM)));
	}

	rep(t,0,20000) {
		int n = (int)rnd(0, 40), mode = (int)rnd(0, 5);
		switch (rnd(0, 4)) {
			case 0: check(gen<int>(n, -5, 5, mode)); break;
			case 1: check(gen<int>(n, -1000000000, 1000000000, mode)); break;
			case 2: check(gen<int>(n, IMIN, IMAX, mode)); break;
			case 3: check(gen<ll>(n, -LIM, LIM, mode)); break;
			default: {
				// int and ll inputs must give identical edges
				auto a = gen<int>(n, IMIN, IMAX, mode);
				vector<Point<ll>> b;
				for (auto p : a) b.push_back(Point<ll>(p.x, p.y));
				assert(manhattanMST(a) == manhattanMST(b));
				check(b);
			}
		}
	}
	rep(t,0,20) check(gen<ll>(1000, -LIM, LIM, t % 6));

	{ // many tiny testcases: cost must scale with n only
		auto start = chrono::steady_clock::now();
		ll sum = 0;
		rep(t,0,200000) {
			vector<Point<int>> a{Point<int>(t, -t), Point<int>(-t, 3), Point<int>(7, t)};
			for (auto e : manhattanMST(a)) sum += e[0];
		}
		double s = chrono::duration<double>(chrono::steady_clock::now() - start).count();
		assert(sum > 0);
		cerr << "200000 tiny cases: " << s << " s" << endl;
		assert(s < 20); // generous: ~0.1 s at -O2, sanitizer builds are slower
	}
	cout << "Tests passed!" << endl;
}
