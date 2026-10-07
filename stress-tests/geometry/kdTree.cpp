#include "../utilities/template.h"

#include "../../content/geometry/kdTree.h"

mt19937_64 rng(12345);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }

void check(const vector<P>& ps, const vector<P>& qs) {
	KDTree tree(ps);
	for (P q : qs) {
		T best = INF;
		for (P p : ps) best = min(best, (p - q).dist2());
		auto res = tree.nearest(q);
		assert(res.first == best);
		assert((res.second - q).dist2() == best);
		assert(count(all(ps), res.second));
	}
}

int main() {
	// tiny grids: lots of duplicates, collinear points and ties
	rep(it,0,60000) {
		int n = (int)rnd(1, 12), r = (int)rnd(0, 4);
		int mode = (int)rnd(0, 3);
		vector<P> ps, qs;
		rep(i,0,n) {
			ll x = rnd(-r, r), y = rnd(-r, r);
			if (mode == 1) y = 0; // horizontal line
			if (mode == 2) x = 0; // vertical line
			if (mode == 3) y = x; // diagonal
			ps.push_back(P(x, y));
		}
		rep(i,0,6) qs.push_back(P(rnd(-r-2, r+2), rnd(-r-2, r+2)));
		for (P p : ps) qs.push_back(p);
		check(ps, qs);
	}
	// medium random, several coordinate ranges up to 1e9
	// (squared distances up to 8e18 still fit in a long long)
	for (ll r : {10LL, 1000LL, 1000000LL, 1000000000LL}) rep(it,0,300) {
		int n = (int)rnd(1, 300);
		vector<P> ps, qs;
		rep(i,0,n) ps.push_back(P(rnd(-r, r), rnd(-r, r)));
		rep(i,0,50) qs.push_back(P(rnd(-r, r), rnd(-r, r)));
		check(ps, qs);
	}
	// extreme corners
	{
		ll r = 1000000000;
		vector<P> ps, qs;
		for (ll x : {-r, r}) for (ll y : {-r, r}) {
			ps.push_back(P(x, y)); qs.push_back(P(x, y));
			rep(i,0,20) ps.push_back(P(x - sgn(x) * rnd(0, 5), y - sgn(y) * rnd(0, 5)));
		}
		qs.push_back(P(0, 0));
		check(ps, qs);
	}
	// all points equal; points on a circle queried from the center
	{
		vector<P> ps(2000, P(7, -3)), qs;
		rep(i,0,20) qs.push_back(P(rnd(-10, 10), rnd(-10, 10)));
		check(ps, qs);
		ps.clear();
		rep(i,0,2000) {
			double a = i * 2 * acos(-1) / 2000;
			ps.push_back(P(llround(1e6 * cos(a)), llround(1e6 * sin(a))));
		}
		qs.push_back(P(0, 0));
		check(ps, qs);
	}
	// large
	{
		int n = 200000;
		vector<P> ps;
		rep(i,0,n) ps.push_back(P(rnd(-1000000, 1000000), rnd(-1000000, 1000000)));
		KDTree tree(ps);
		rep(it,0,200) {
			P q(rnd(-1100000, 1100000), rnd(-1100000, 1100000));
			T best = INF;
			for (P p : ps) best = min(best, (p - q).dist2());
			assert(tree.nearest(q).first == best);
		}
	}
	// Performance regression: collinear points on the diagonal y = x with
	// random queries. A search that does not carry the best distance found so
	// far into the far subtree visits ~n/9 nodes per query here (5 s for this
	// block); a standard one visits O(sqrt n) (0.2 s).
	{
		int n = 1000000; ll r = 1000000000;
		vector<ll> xs(n);
		for (ll& x : xs) x = rnd(-r, r);
		sort(all(xs));
		vector<P> ps;
		for (ll x : xs) ps.push_back(P(x, x));
		shuffle(all(ps), rng);
		KDTree tree(ps);
		vector<P> qs;
		rep(i,0,3000) qs.push_back(P(rnd(-r, r), rnd(-r, r)));
		vector<T> got;
		auto start = chrono::steady_clock::now();
		for (P q : qs) got.push_back(tree.nearest(q).first);
		double secs = chrono::duration<double>(chrono::steady_clock::now() - start).count();
		rep(i,0,sz(qs)) { // dist2 is convex in x: only look around (qx+qy)/2
			P q = qs[i];
			int j = int(lower_bound(all(xs), (q.x + q.y) / 2) - xs.begin());
			T best = INF;
			rep(k,max(0,j-2),min(n,j+3)) best = min(best, (P(xs[k], xs[k]) - q).dist2());
			assert(got[i] == best);
		}
		if (secs > 1.0) {
			cerr << "kdTree: 3000 queries on 1e6 collinear points took " << secs << " s" << endl;
			return 1;
		}
	}
	cout<<"Tests passed!"<<endl;
}
