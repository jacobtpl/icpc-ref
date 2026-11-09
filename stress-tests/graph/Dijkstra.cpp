#include "../utilities/template.h"

#include "../../content/graph/Dijkstra.h"
#include "../../content/graph/BellmanFord.h"

int main() {
	rep(it,0,3000) {
		int n = rand() % 8 + 1, m = rand() % 15;
		vector<vector<pair<int, ll>>> g(n);
		vector<Ed> eds;
		rep(i,0,m) {
			int a = rand() % n, b = rand() % n, w = rand() % 10;
			g[a].emplace_back(b, w);
			eds.push_back({a, b, w});
		}
		int s = rand() % n;
		vector<Node> nodes(n);
		bellmanFord(nodes, eds, s);

		// plain wrapper
		auto d = dijkstra(g, s);
		rep(i,0,n) assert(d[i] == nodes[i].dist);

		// generic version with parent reconstruction, D = pair (dist, hops)
		using D = pair<ll, int>;
		auto r = dijkstra<D>(n, {{s, D(0, 0)}}, less<D>(),
			[&](int u, D du, auto push) {
				for (auto [v, w] : g[u]) push(v, D(du.first + w, du.second + 1));
			});
		rep(i,0,n) {
			assert(r.seen[i] == (nodes[i].dist != inf));
			if (!r.seen[i]) continue;
			assert(r.dist[i].first == nodes[i].dist);
			// walk parents back to s; check path is consistent
			int v = i, hops = 0; ll len = 0;
			while (r.par[v] != -1) {
				int p = r.par[v];
				ll best = LLONG_MAX;
				for (auto [x, w] : g[p]) if (x == v) best = min(best, w);
				assert(best != LLONG_MAX);
				len += best; hops++; v = p;
			}
			assert(v == s && len == r.dist[i].first && hops == r.dist[i].second);
		}

		// multi-source, implicit edges (generated on the fly)
		int s2 = rand() % n;
		auto r2 = dijkstra<ll>(n, {{s, 0}, {s2, 3}}, less<ll>(),
			[&](int u, ll du, auto push) {
				for (auto [v, w] : g[u]) push(v, du + w);
			});
		vector<Node> n1(n), n2(n);
		bellmanFord(n1, eds, s); bellmanFord(n2, eds, s2);
		rep(i,0,n) {
			ll e = min(n1[i].dist, n2[i].dist == inf ? inf : n2[i].dist + 3);
			assert((e != inf) == (bool)r2.seen[i]);
			if (r2.seen[i]) assert(r2.dist[i] == e);
		}
	}
	// large test: path graph, 2e5 nodes
	{
		int n = 200000;
		vector<vector<pair<int, ll>>> g(n);
		rep(i,0,n-1) g[i].emplace_back(i + 1, 1000000000);
		auto d = dijkstra(g, 0);
		assert(d[n-1] == (ll)(n - 1) * 1000000000);
	}
	cout << "Tests passed!" << endl;
}
