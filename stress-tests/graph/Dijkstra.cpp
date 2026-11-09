#include "../utilities/template.h"

#include "../../content/graph/Dijkstra.h"
#include "../../content/graph/BellmanFord.h"

// custom distance: lexicographic (length, #edges)
struct HopsDist {
	typedef pair<ll, int> D;
	D zero() { return {0, 0}; }
	D add(D d, ll w) { return {d.first + w, d.second + 1}; }
	bool less(D a, D b) { return a < b; }
};

int main() {
	rep(it,0,3000) {
		int n = rand() % 8 + 1, m = rand() % 15;
		vector<vector<pair<int, ll>>> g(n);
		vector<vector<pii>> gi(n); // int weights -> default LLDist still works
		vector<Ed> eds;
		rep(i,0,m) {
			int a = rand() % n, b = rand() % n, w = rand() % 10;
			g[a].emplace_back(b, w);
			gi[a].emplace_back(b, w);
			eds.push_back({a, b, w});
		}
		int s = rand() % n;
		vector<Node> nodes(n);
		bellmanFord(nodes, eds, s);

		// default distance (ll)
		auto r0 = dijkstra(g, {s});
		auto r1 = dijkstra(gi, {s});
		rep(i,0,n) {
			assert((bool)r0.seen[i] == (nodes[i].dist != inf));
			assert(r1.seen[i] == r0.seen[i]);
			if (r0.seen[i]) assert(r0.dist[i] == nodes[i].dist && r1.dist[i] == nodes[i].dist);
		}

		// custom distance with parent reconstruction
		auto r = dijkstra(g, {s}, HopsDist{});
		rep(i,0,n) {
			assert(r.seen[i] == r0.seen[i]);
			if (!r.seen[i]) continue;
			assert(r.dist[i].first == nodes[i].dist);
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

		// multi-source
		int s2 = rand() % n;
		auto r2 = dijkstra(g, {s, s2});
		vector<Node> n1(n), n2(n);
		bellmanFord(n1, eds, s); bellmanFord(n2, eds, s2);
		rep(i,0,n) {
			ll e = min(n1[i].dist, n2[i].dist);
			assert((e != inf) == (bool)r2.seen[i]);
			if (r2.seen[i]) assert(r2.dist[i] == e);
		}

		// implicit graph (edges generated on the fly)
		auto r3 = dijkstraGen(n, {s}, LLDist{}, [&](int u, auto f) {
			for (auto [v, w] : g[u]) f(v, w);
		});
		rep(i,0,n) {
			assert(r3.seen[i] == r0.seen[i]);
			if (r0.seen[i]) assert(r3.dist[i] == r0.dist[i]);
		}
	}
	// large test: path graph, 2e5 nodes
	{
		int n = 200000;
		vector<vector<pair<int, ll>>> g(n);
		rep(i,0,n-1) g[i].emplace_back(i + 1, 1000000000);
		auto r = dijkstra(g, {0});
		assert(r.dist[n-1] == (ll)(n - 1) * 1000000000);
	}
	cout << "Tests passed!" << endl;
}
