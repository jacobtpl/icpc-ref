#include "../utilities/template.h"
#include "../utilities/random.h"

#include "../../content/graph/DominatorTree.h"

// Brute force: v is reachable from r; u strictly dominates v iff v is
// unreachable from r once u is removed. idom = the strict dominator
// with the most dominators (the deepest one).
vi reach(vector<vi>& g, int r, int del) {
	int n = sz(g);
	vi seen(n);
	if (r == del) return seen;
	vi st{r}; seen[r] = 1;
	while (!st.empty()) {
		int v = st.back(); st.pop_back();
		for (int u : g[v]) if (u != del && !seen[u]) seen[u] = 1, st.push_back(u);
	}
	return seen;
}
vi brute(vector<vi>& g, int r) {
	int n = sz(g);
	vi base = reach(g, r, -1), res(n, -1);
	vector<vi> doms(n);
	rep(u,0,n) if (base[u]) {
		vi s = reach(g, r, u);
		rep(v,0,n) if (base[v] && !s[v] && v != u) doms[v].push_back(u);
	}
	res[r] = r;
	rep(v,0,n) if (base[v] && v != r) {
		int best = -1;
		for (int u : doms[v])
			if (best == -1 || sz(doms[u]) > sz(doms[best])) best = u;
		assert(best != -1);
		res[v] = best;
	}
	return res;
}

void check(vector<vi>& g, int r) {
	vi got = domTree(g, r), want = brute(g, r);
	if (got != want) {
		cerr << "n=" << sz(g) << " r=" << r << endl;
		rep(i,0,sz(g)) for (int j : g[i]) cerr << i << "->" << j << endl;
		for (int x : got) cerr << x << ' ';
		cerr << endl;
		for (int x : want) cerr << x << ' ';
		cerr << endl;
		assert(0);
	}
}

int main() {
	// exhaustive: all digraphs (with self-loops) on n <= 3, all roots
	rep(n,1,4) rep(bits,0,1 << (n*n)) {
		vector<vi> g(n);
		rep(i,0,n) rep(j,0,n) if (bits >> (i*n+j) & 1) g[i].push_back(j);
		rep(r,0,n) check(g, r);
	}
	// random small multigraphs
	rep(it,0,100000) {
		int n = rand() % 10 + 1, m = rand() % 25;
		vector<vi> g(n);
		rep(i,0,m) g[rand() % n].push_back(rand() % n);
		check(g, rand() % n);
	}
	// random medium graphs, several shapes
	rep(it,0,1000) {
		int n = rand() % 150 + 1, type = rand() % 3;
		vector<vi> g(n);
		if (type == 0) {
			int m = rand() % (3 * n + 1);
			rep(i,0,m) g[rand() % n].push_back(rand() % n);
		} else if (type == 1) { // DAG-ish: edges mostly forward
			int m = rand() % (3 * n + 1);
			rep(i,0,m) {
				int a = rand() % n, b = rand() % n;
				if (a > b && rand() % 4) swap(a, b);
				g[a].push_back(b);
			}
		} else { // path with extra back/forward edges
			rep(i,1,n) g[i-1].push_back(i);
			int m = rand() % (n + 1);
			rep(i,0,m) g[rand() % n].push_back(rand() % n);
		}
		check(g, rand() % n);
	}
	// edge cases
	{
		vector<vi> g(1);
		assert(domTree(g, 0) == vi{0});
		g[0].push_back(0); g[0].push_back(0);
		assert(domTree(g, 0) == vi{0});
		vector<vi> h(3); // isolated root, others unreachable
		h[1].push_back(2); h[2].push_back(1);
		assert(domTree(h, 0) == (vi{0, -1, -1}));
	}
	// long path with back edges: no recursion, O((V+E) log V)
	{
		int n = 200000;
		vector<vi> g(n);
		rep(i,1,n) g[i-1].push_back(i);
		rep(i,0,n) g[i].push_back(rand() % (i + 1));
		vi d = domTree(g, 0);
		rep(i,1,n) assert(d[i] == i - 1);
		assert(d[0] == 0);
		// fan from last vertex: everything reachable, idom of fan = n-1
		vector<vi> h(2 * n);
		rep(i,1,n) h[i-1].push_back(i);
		rep(i,n,2*n) h[n-1].push_back(i), h[i].push_back(rand() % (2 * n));
		d = domTree(h, 0);
		rep(i,n,2*n) assert(d[i] == n - 1);
	}
	cout << "Tests passed!" << endl;
}
