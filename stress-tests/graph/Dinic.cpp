#include "../utilities/template.h"

#include "../../content/graph/Dinic.h"
// Pasted into the same solution, the two headers must not clash: both used
// to define a global `struct Edge`.
#include "../../content/graph/DirectedMST.h"

mt19937_64 rng(4242);
ll rnd(ll lo, ll hi) { return (ll)(rng() % (unsigned long long)(hi - lo + 1)) + lo; }

// Oracle: Ford-Fulkerson with DFS on an adjacency matrix (parallel edges merged).
struct Brute {
	int n;
	vector<vector<ll>> c;
	vi vis;
	Brute(int n) : n(n), c(n, vector<ll>(n)) {}
	ll dfs(int u, int t, ll f) {
		if (u == t) return f;
		vis[u] = 1;
		rep(v,0,n) if (!vis[v] && c[u][v] > 0)
			if (ll p = dfs(v, t, min(f, c[u][v]))) {
				c[u][v] -= p, c[v][u] += p;
				return p;
			}
		return 0;
	}
	ll calc(int s, int t) {
		ll flow = 0;
		for (;;) {
			vis.assign(n, 0);
			ll p = dfs(s, t, LLONG_MAX);
			if (!p) return flow;
			flow += p;
		}
	}
};

struct E3 { int a, b; ll c; };

void check(int n, int s, int t, const vector<E3>& eds) {
	Dinic d(n);
	Brute b(n);
	int real = 0;
	for (auto& e : eds) {
		d.AddEdge(e.a, e.b, e.c);
		if (e.a != e.b) b.c[e.a][e.b] += e.c, real++;
	}
	assert(sz(d.E) == 2 * real);
	ll flow = d.MaxFlow(s, t), expected = b.calc(s, t);
	assert(flow == expected);
	// capacity constraints, antisymmetry and flow conservation
	vector<ll> ex(n);
	rep(i,0,sz(d.E)) {
		auto& e = d.E[i];
		assert(e.flow <= e.cap);
		assert(e.flow == -d.E[i ^ 1].flow);
		if (i % 2 == 0) {
			assert(e.flow >= 0);
			ex[e.u] -= e.flow, ex[e.v] += e.flow;
		}
	}
	rep(i,0,n) if (i != s && i != t) assert(ex[i] == 0);
	assert(ex[t] == flow && ex[s] == -flow);
	// the last BFS labels define a minimum cut
	ll cut = 0;
	assert(d.d[s] == 0 && d.d[t] == n + 1);
	rep(i,0,sz(d.E)) if (i % 2 == 0) {
		auto& e = d.E[i];
		if (d.d[e.u] <= n && d.d[e.v] > n) cut += e.cap;
	}
	assert(cut == flow);
}

int main() {
	// tiny graphs, including self-loops, multi-edges and zero capacities
	rep(it,0,300000) {
		int n = (int)rnd(2, 8);
		int s = (int)rnd(0, n-1), t = (int)rnd(0, n-2);
		if (t >= s) t++;
		int m = (int)rnd(0, it % 2 ? 30 : 8);
		ll maxc = it % 5 == 0 ? 1 : it % 5 == 1 ? 1000000000 : 5;
		vector<E3> eds;
		rep(i,0,m) eds.push_back({(int)rnd(0, n-1), (int)rnd(0, n-1), rnd(0, maxc)});
		check(n, s, t, eds);
	}
	// medium graphs
	rep(it,0,300) {
		int n = (int)rnd(20, 60);
		int m = (int)rnd(0, n * 6);
		vector<E3> eds;
		rep(i,0,m) eds.push_back({(int)rnd(0, n-1), (int)rnd(0, n-1), rnd(0, 100)});
		check(n, 0, n - 1, eds);
	}
	// layered graphs (many blocking-flow phases, long paths)
	rep(it,0,300) {
		int L = (int)rnd(2, 8), W = (int)rnd(1, 5), n = L * W + 2;
		vector<E3> eds;
		rep(j,0,W) eds.push_back({0, 1 + j, rnd(0, 20)});
		rep(l,0,L-1) rep(i,0,W) rep(j,0,W) if (rnd(0, 2))
			eds.push_back({1 + l * W + i, 1 + (l + 1) * W + j, rnd(0, 20)});
		rep(j,0,W) eds.push_back({1 + (L - 1) * W + j, n - 1, rnd(0, 20)});
		rep(k,0,(int)rnd(0, 5)) eds.push_back({(int)rnd(1, n-2), (int)rnd(1, n-2), rnd(0, 20)});
		shuffle(all(eds), rng);
		check(n, 0, n - 1, eds);
	}
	// huge capacities: 20 edges of up to 4e17 keep every sum below LLONG_MAX
	rep(it,0,20000) {
		int n = (int)rnd(2, 7);
		int m = (int)rnd(0, 20);
		vector<E3> eds;
		const ll big = 400000000000000000LL;
		rep(i,0,m) eds.push_back({(int)rnd(0, n-1), (int)rnd(0, n-1), rnd(big - 5, big)});
		check(n, 0, n - 1, eds);
	}
	// a single edge of capacity LLONG_MAX
	check(2, 0, 1, {{0, 1, LLONG_MAX}});
	check(3, 0, 2, {{0, 1, LLONG_MAX}, {1, 2, LLONG_MAX}});
	// disconnected sink, no edges
	check(2, 0, 1, {});
	check(5, 0, 4, {{0, 1, 3}, {1, 2, 3}, {3, 4, 3}});
	// long path: recursion depth n
	{
		int n = 30000;
		vector<E3> eds;
		rep(i,0,n-1) eds.push_back({i, i + 1, 7 + i % 3});
		Dinic d(n);
		for (auto& e : eds) d.AddEdge(e.a, e.b, e.c);
		assert(d.MaxFlow(0, n - 1) == 7);
	}
	cout << "Tests passed!" << endl;
}
