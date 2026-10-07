#include "../utilities/template.h"

#include "../../content/graph/PushRelabel.h"

// Oracle: Ford-Fulkerson with DFS on an adjacency matrix.
struct FF {
	int n; vector<vector<ll>> c; vi vis;
	FF(int n) : n(n), c(n, vector<ll>(n)) {}
	ll dfs(int u, int t, ll f) {
		if (u == t) return f;
		vis[u] = 1;
		rep(v,0,n) if (!vis[v] && c[u][v] > 0)
			if (ll d = dfs(v, t, min(f, c[u][v]))) return c[u][v] -= d, c[v][u] += d, d;
		return 0;
	}
	ll calc(int s, int t) {
		ll r = 0;
		for (ll d;; r += d) { vis.assign(n, 0); if (!(d = dfs(s, t, LLONG_MAX))) return r; }
	}
};

mt19937 rng(777);
int rnd(int lo, int hi) { return lo + (int)(rng() % (unsigned)(hi - lo + 1)); }

struct Ed { int a, b; ll c, rc; };
void check(int n, int s, int t, const vector<Ed>& es) {
	PushRelabel pr(n); FF ff(n);
	vector<vector<ll>> cap(n, vector<ll>(n));
	for (auto& e : es) {
		pr.addEdge(e.a, e.b, e.c, e.rc);
		if (e.a != e.b) ff.c[e.a][e.b] += e.c, ff.c[e.b][e.a] += e.rc, cap[e.a][e.b] += e.c, cap[e.b][e.a] += e.rc;
	}
	ll flow = pr.calc(s, t);
	assert(flow == ff.calc(s, t));
	// "look at positive values only" gives a valid max flow
	vector<ll> bal(n);
	rep(i,0,n) for (auto& e : pr.g[i]) {
		assert(e.c >= 0);
		assert(e.f == -pr.g[e.dest][e.back].f);
		if (e.f > 0) bal[i] += e.f, bal[e.dest] -= e.f;
	}
	rep(i,0,n) assert(bal[i] == (i == s ? flow : i == t ? -flow : 0));
	// leftOfMinCut is a minimum cut
	assert(pr.leftOfMinCut(s) && !pr.leftOfMinCut(t));
	ll cut = 0;
	rep(i,0,n) rep(j,0,n) if (pr.leftOfMinCut(i) && !pr.leftOfMinCut(j)) cut += cap[i][j];
	assert(cut == flow);
}

int main() {
	check(2, 0, 1, {});
	check(2, 0, 1, {{0, 1, 5, 0}});
	check(2, 1, 0, {{0, 1, 5, 0}});
	check(2, 0, 1, {{0, 1, 5, 0}, {0, 1, 7, 2}, {0, 0, 3, 3}, {1, 0, 4, 0}});
	check(3, 0, 2, {{0, 1, (ll)4e18, 0}, {1, 2, (ll)3e18, 0}}); // near ll limit
	check(4, 0, 3, {{0, 1, (ll)2e18, 0}, {0, 2, (ll)2e18, 0}, {1, 3, (ll)2e18, 0}, {2, 3, (ll)2e18, 0}});
	rep(it,0,300000) {
		int n = rnd(2, 8), m = rnd(0, 25), s = rnd(0, n - 1), t = rnd(0, n - 2), mx = it % 3 ? 4 : 1000000000;
		if (t >= s) t++;
		vector<Ed> es;
		rep(i,0,m) es.push_back({rnd(0, n - 1), rnd(0, n - 1), rnd(0, mx), rnd(0, 3) ? 0 : rnd(0, mx)});
		check(n, s, t, es);
	}
	rep(it,0,3000) {
		int n = rnd(2, 60), m = rnd(0, 6 * n), s = rnd(0, n - 1), t = rnd(0, n - 2);
		if (t >= s) t++;
		vector<Ed> es;
		rep(i,0,m) es.push_back({rnd(0, n - 1), rnd(0, n - 1), rnd(0, it % 2 ? 3 : 1000), rnd(0, 1) ? 0 : rnd(0, 5)});
		check(n, s, t, es);
	}
	// layered / path-like graphs
	rep(it,0,2000) {
		int L = rnd(1, 8), W = rnd(1, 5), n = L * W + 2;
		vector<Ed> es;
		rep(j,0,W) es.push_back({n - 2, j, rnd(0, 9), 0}), es.push_back({(L - 1) * W + j, n - 1, rnd(0, 9), 0});
		rep(l,0,L-1) rep(i,0,W) rep(j,0,W) if (rnd(0, 1)) es.push_back({l * W + i, (l + 1) * W + j, rnd(0, 9), 0});
		check(n, n - 2, n - 1, es);
	}
	cout << "Tests passed!" << endl;
}
