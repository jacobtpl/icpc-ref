#include "../utilities/template.h"

#include "../../content/graph/BiconnectedComponents.h"

mt19937 rng(31337);
int ri(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }

// Oracle: subdivide every edge; two edges are in the same block iff their
// midpoints are connected and stay connected after deleting any one vertex.
// Returns the classes with >= 2 edges (singletons are bridges), sorted.
vector<vi> oracle(int n, const vector<pii>& es) {
	int m = sz(es);
	vector<vi> same(m, vi(m, 1));
	rep(del,-1,n) {
		vi comp(n + m, -1);
		int c = 0;
		rep(st,0,n+m) if (comp[st] == -1 && st != del) {
			vi q = {st};
			comp[st] = c;
			rep(i,0,sz(q)) {
				int v = q[i];
				auto go = [&](int y) {
					if (y != del && comp[y] == -1) comp[y] = c, q.push_back(y);
				};
				if (v >= n) go(es[v-n].first), go(es[v-n].second);
				else rep(e,0,m) if (es[e].first == v || es[e].second == v) go(n + e);
			}
			c++;
		}
		rep(a,0,m) rep(b,0,m) if (comp[n+a] != comp[n+b]) same[a][b] = 0;
	}
	vector<vi> res;
	rep(a,0,m) {
		vi cl;
		rep(b,0,m) if (same[a][b]) cl.push_back(b);
		if (cl[0] == a && sz(cl) > 1) res.push_back(cl);
	}
	return res;
}

vector<vi> run(int n, const vector<pii>& es) {
	ed.assign(n, {});
	rep(i,0,sz(es)) {
		ed[es[i].first].emplace_back(es[i].second, i);
		ed[es[i].second].emplace_back(es[i].first, i);
	}
	for (auto& v : ed) shuffle(all(v), rng);
	vector<vi> res;
	bicomps([&](const vi& edgelist) {
		res.push_back(edgelist);
		sort(all(res.back()));
	});
	assert(st.empty());
	sort(all(res));
	return res;
}

int main() {
	rep(it,0,150000) {
		int n = ri(1, 8), m = n == 1 ? 0 : ri(0, it % 4 ? 2*n : 3);
		vector<pii> es;
		rep(i,0,m) { // no self-loops; multi-edges allowed
			int a = ri(0, n-1), b = ri(0, n-2);
			if (b >= a) b++;
			es.emplace_back(a, b);
		}
		auto got = run(n, es), want = oracle(n, es);
		if (got != want) {
			cerr << "mismatch n=" << n << endl;
			for (auto e : es) cerr << e.first << ' ' << e.second << endl;
			abort();
		}
	}
	// empty graph, no edges, repeated calls on the same globals
	assert(run(0, {}).empty());
	assert(run(5, {}).empty());
	assert(run(2, {{0, 1}}).empty());
	assert(run(2, {{0, 1}, {1, 0}}) == vector<vi>({{0, 1}}));
	// larger structured graphs: cactus-like chains of cycles joined by bridges
	rep(it,0,300) {
		int n = 0;
		vector<pii> es;
		vector<vi> want;
		int blocks = ri(1, 30), last = 0;
		n = 1;
		rep(b,0,blocks) {
			int at = ri(0, n-1), len = ri(1, 6);
			(void)last;
			if (len == 1) { es.emplace_back(at, n++); continue; } // bridge
			vi cyc;
			int prev = at;
			rep(i,1,len) cyc.push_back(sz(es)), es.emplace_back(prev, n), prev = n++;
			cyc.push_back(sz(es)), es.emplace_back(prev, at);
			want.push_back(cyc);
		}
		vi perm(n);
		iota(all(perm), 0);
		shuffle(all(perm), rng);
		for (auto& e : es) e = {perm[e.first], perm[e.second]};
		sort(all(want));
		assert(run(n, es) == want);
	}
	cout << "Tests passed!" << endl;
}
