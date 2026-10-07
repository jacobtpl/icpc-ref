#include "../utilities/template.h"

#include "../../content/graph/BellmanFord.h"

mt19937 rng(12345);
int ri(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }

struct E { int a, b, w; };
// Oracle: textbook Bellman-Ford with n full rounds, then propagate -inf
// from every node that can still be relaxed.
vector<ll> oracle(int n, const vector<E>& es, int s) {
	vector<ll> d(n, inf);
	d[s] = 0;
	rep(it,0,n) for (auto& e : es) if (d[e.a] != inf)
		d[e.b] = min(d[e.b], d[e.a] + e.w);
	vector<bool> neg(n);
	for (auto& e : es) if (d[e.a] != inf && d[e.a] + e.w < d[e.b]) neg[e.b] = 1;
	rep(it,0,n) for (auto& e : es) if (neg[e.a]) neg[e.b] = 1;
	rep(i,0,n) if (neg[i]) d[i] = -inf;
	return d;
}

void check(int n, const vector<E>& es, int s) {
	vector<Node> nodes(n);
	vector<Ed> eds;
	for (auto& e : es) eds.push_back({e.a, e.b, e.w});
	bellmanFord(nodes, eds, s);
	vector<ll> d = oracle(n, es, s);
	rep(i,0,n) {
		if (nodes[i].dist != d[i]) {
			cerr << "mismatch n=" << n << " s=" << s << " node " << i << " got "
				<< nodes[i].dist << " want " << d[i] << endl;
			for (auto& e : es) cerr << e.a << ' ' << e.b << ' ' << e.w << endl;
			abort();
		}
		// prev must be a tight incoming edge for finitely-reachable nodes
		if (i != s && d[i] != inf && d[i] != -inf) {
			int p = nodes[i].prev;
			assert(p != -1 && d[p] != inf && d[p] != -inf);
			bool ok = 0;
			for (auto& e : es) if (e.a == p && e.b == i && d[p] + e.w == d[i]) ok = 1;
			assert(ok);
		}
	}
	if (d[s] == 0) assert(nodes[s].prev == -1);
}

int main() {
	// tiny graphs, small weights: lots of negative cycles, self-loops, multi-edges
	rep(it,0,300000) {
		int n = ri(1, 7), m = ri(0, 12), lo = ri(-5, 0), hi = ri(0, 8);
		vector<E> es;
		rep(i,0,m) es.push_back({ri(0, n-1), ri(0, n-1), ri(lo, hi)});
		check(n, es, ri(0, n-1));
	}
	// mostly positive weights, few negatives: long shortest paths, no neg cycle
	rep(it,0,20000) {
		int n = ri(1, 30), m = ri(0, 3*n);
		vi pot(n);
		for (int& x : pot) x = ri(-50, 50);
		vector<E> es;
		rep(i,0,m) {
			int a = ri(0, n-1), b = ri(0, n-1);
			es.push_back({a, b, ri(0, 20) + pot[a] - pot[b]});
		}
		check(n, es, ri(0, n-1));
	}
	// adversarial orders for the n/2+2 round limit: paths / zig-zags under
	// random relabelling, optionally closed into a negative cycle at the end
	rep(it,0,20000) {
		int n = ri(2, 40);
		vi perm(n);
		iota(all(perm), 0);
		int kind = ri(0, 3);
		if (kind == 1) reverse(all(perm));
		if (kind == 2) shuffle(all(perm), rng);
		if (kind == 3) { // zig-zag: 0, n-1, 1, n-2, ...
			int l = 0, r = n-1;
			rep(i,0,n) perm[i] = i % 2 ? r-- : l++;
		}
		vector<E> es;
		rep(i,0,n-1) es.push_back({perm[i], perm[i+1], ri(-3, 3)});
		int extra = ri(0, 2);
		if (extra == 1) es.push_back({perm[n-1], perm[ri(0, n-1)], -1000});
		if (extra == 2) es.push_back({perm[n-1], perm[n-1], -1});
		shuffle(all(es), rng);
		check(n, es, perm[ri(0, 1) ? 0 : ri(0, n-1)]);
	}
	// weights at the int limits (documented: V^2 max|w| < 2^63)
	rep(it,0,20000) {
		int n = ri(1, 50), m = ri(0, 4*n);
		vector<E> es;
		rep(i,0,m) {
			int w = ri(0, 3) ? INT_MAX - ri(0, 3) : INT_MIN + ri(0, 3);
			if (ri(0, 9) == 0) w = ri(-2, 2);
			es.push_back({ri(0, n-1), ri(0, n-1), w});
		}
		check(n, es, ri(0, n-1));
	}
	{ // n = 1, no edges
		check(1, {}, 0);
		check(1, {{0, 0, 5}}, 0);
		check(1, {{0, 0, -5}}, 0);
	}
	cout << "Tests passed!" << endl;
}
