// Regression tests for kdTree.h: ownership (the old tree leaked every node;
// run under ASan/LSan to check frees, copies and moves), nearest neighbour
// against brute force over the documented coordinate range, duplicates,
// negative coordinates, empty trees and many tiny testcases.
#include "../utilities/template.h"

#include "../../content/geometry/Point.h"
#include "../../content/geometry/kdTree.h"

static_assert(is_copy_constructible<KDTree>::value && is_copy_assignable<KDTree>::value);
static_assert(is_nothrow_move_constructible<KDTree>::value &&
	is_nothrow_move_assignable<KDTree>::value);
static_assert(is_default_constructible<KDTree>::value);

typedef __int128 lll;
mt19937_64 rng(20261012);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }
const ll LIM = (1LL << 30) - 1; // documented bound

lll d2(P a, P b) {
	lll dx = (lll)a.x - b.x, dy = (lll)a.y - b.y;
	return dx * dx + dy * dy;
}

void checkQuery(const KDTree& tree, const vector<P>& ps, P q) {
	pair<T, P> r = tree.nearest(q);
	if (ps.empty()) {
		assert(r.first == INF && r.second == P());
		return;
	}
	lll best = d2(ps[0], q);
	for (P p : ps) best = min(best, d2(p, q));
	assert(r.first == best);
	assert(d2(r.second, q) == best);
	assert(find(all(ps), r.second) != ps.end());
}

vector<P> gen(int n, ll lo, ll hi, int mode) {
	vector<P> ps;
	rep(i,0,n) {
		ll x = rnd(lo, hi), y = rnd(lo, hi);
		if (mode == 1) x = rnd(0, 1) ? lo : hi, y = rnd(0, 1) ? lo : hi; // corners
		if (mode == 2) y = x;              // diagonal
		if (mode == 3) y = lo;             // horizontal line
		if (mode == 4) x = hi;             // vertical line
		if (mode == 5 && i && rnd(0, 2)) ps.push_back(ps[rnd(0, i - 1)]); // duplicates
		else ps.push_back(P(x, y));
	}
	return ps;
}

P genQuery(const vector<P>& ps, ll lo, ll hi) {
	switch (rnd(0, 3)) {
		case 0: return P(rnd(0, 1) ? lo : hi, rnd(0, 1) ? lo : hi);
		case 1: if (!ps.empty()) return ps[rnd(0, sz(ps) - 1)]; [[fallthrough]];
		case 2: if (!ps.empty()) {
			P p = ps[rnd(0, sz(ps) - 1)];
			return P(clamp(p.x + rnd(-2, 2), lo, hi), clamp(p.y + rnd(-2, 2), lo, hi));
		} [[fallthrough]];
		default: return P(rnd(lo, hi), rnd(lo, hi));
	}
}

int main() {
	{ // the review reproducer (leaked 960 bytes in 15 allocations)
		for (int i = 0; i < 3; ++i) {
			KDTree tree({P(0,0), P(2,2), P(4,0)});
			assert(tree.nearest(P(0,0)).first == 0);
			assert(tree.nearest(P(3,3)) == make_pair(2LL, P(2,2)));
		}
	}
	{ // maximal squared distance in the documented range: 8 * (2^30-1)^2
		KDTree tree({P(-LIM, -LIM)});
		auto r = tree.nearest(P(LIM, LIM));
		assert(r.first == 8 * LIM * LIM && r.second == P(-LIM, -LIM));
		KDTree two({P(-LIM, -LIM), P(-LIM, LIM)});
		assert(two.nearest(P(LIM, -LIM)).first == 4 * LIM * LIM);
		assert(two.nearest(P(LIM, 1)).second == P(-LIM, LIM));
	}
	{ // empty trees
		KDTree a, b((vector<P>()));
		assert(a.nearest(P(1, 2)) == make_pair(INF, P()));
		assert(b.nearest(P(-1, -2)).first == INF);
		KDTree c = a; c = b; b = move(c);
		assert(b.nearest(P()).first == INF);
	}
	{ // all points equal
		vector<P> ps(9, P(-7, 3));
		KDTree tree(ps);
		assert(tree.nearest(P(-7, 3)) == make_pair(0LL, P(-7, 3)));
		assert(tree.nearest(P(-4, 7)) == make_pair(25LL, P(-7, 3)));
	}

	// brute force comparison
	rep(t,0,30000) {
		int n = (int)rnd(0, 30), mode = (int)rnd(0, 5);
		ll hi = rnd(0, 2) == 0 ? LIM : rnd(0, 1) ? 3 : 1000, lo = -hi;
		if (rnd(0, 9) == 0) lo = hi - 5; // far from origin
		vector<P> ps = gen(n, lo, hi, mode);
		KDTree tree(ps);
		assert(sz(tree.t) == max(2 * n - 1, 0));
		rep(q,0,20) checkQuery(tree, ps, genQuery(ps, -hi, hi));
	}
	rep(t,0,30) {
		ll hi = t % 2 ? LIM : 100;
		vector<P> ps = gen(2000, -hi, hi, t % 6);
		KDTree tree(ps);
		rep(q,0,200) checkQuery(tree, ps, genQuery(ps, -hi, hi));
	}

	// copies, moves and assignments keep answering correctly and never
	// share or leak storage (checked by ASan/LSan)
	rep(t,0,2000) {
		vector<P> pa = gen((int)rnd(0, 20), -LIM, LIM, (int)rnd(0, 5));
		vector<P> pb = gen((int)rnd(0, 20), -50, 50, (int)rnd(0, 5));
		auto ok = [&](const KDTree& tr, const vector<P>& ps) {
			rep(q,0,5) checkQuery(tr, ps, genQuery(ps, -LIM, LIM));
		};
		auto* a = new KDTree(pa);
		KDTree copy(*a);          // copy construct
		KDTree b(pb);
		KDTree assigned(pb);
		assigned = *a;            // copy assign over a live tree
		delete a;                 // copies must survive the original
		ok(copy, pa); ok(assigned, pa);
		KDTree moved(std::move(copy)); // move construct
		ok(moved, pa);
		copy = b;                 // moved-from tree can be reassigned
		ok(copy, pb);
		b = std::move(moved);     // move assign over a live tree
		ok(b, pa);
		KDTree& self = b;
		b = self;                 // self assignment
		ok(b, pa);
		swap(b, copy);
		ok(b, pb); ok(copy, pa);
		vector<KDTree> v;         // reallocation moves/copies trees
		rep(i,0,9) v.push_back(i % 2 ? KDTree(pa) : KDTree(pb));
		v.erase(v.begin() + 3);
		rep(i,0,sz(v)) ok(v[i], (i < 3 ? i : i + 1) % 2 ? pa : pb);
	}

	{ // many tiny testcases: cost and memory must scale with n only
		auto start = chrono::steady_clock::now();
		ll sum = 0;
		rep(t,0,200000) {
			KDTree tree({P(t, -t), P(-t, 3), P(7, t)});
			sum += tree.nearest(P(1, 1)).first;
		}
		double s = chrono::duration<double>(chrono::steady_clock::now() - start).count();
		assert(sum > 0);
		cerr << "200000 tiny trees: " << s << " s" << endl;
		assert(s < 20); // generous: ~0.05 s at -O2, sanitizer builds are slower
	}
	cout << "Tests passed!" << endl;
}
