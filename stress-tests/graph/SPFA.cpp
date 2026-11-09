#include "../utilities/template.h"

#include "../../content/graph/SPFA.h"
#include "../../content/graph/BellmanFord.h"

int main() {
	int withCycle = 0;
	rep(it,0,20000) {
		int n = rand() % 7 + 1, m = rand() % 14;
		vector<vector<pair<int, ll>>> g(n);
		vector<Ed> eds;
		rep(i,0,m) {
			int a = rand() % n, b = rand() % n, w = rand() % 21 - 10;
			g[a].emplace_back(b, w);
			eds.push_back({a, b, w});
		}
		int s = rand() % n;
		vector<Node> nodes(n);
		bellmanFord(nodes, eds, s);
		bool neg = false;
		rep(i,0,n) neg |= nodes[i].dist == -inf;

		auto r = spfa(g, {s});
		if (neg) {
			withCycle++;
			// must report a cycle: consecutive vertices joined by edges, total < 0,
			// and reachable from s
			auto& c = r.negCycle;
			assert(!c.empty());
			ll tot = 0;
			rep(i,0,sz(c)) {
				int a = c[i], b = c[(i + 1) % sz(c)];
				ll best = LLONG_MAX;
				for (auto [v, w] : g[a]) if (v == b) best = min(best, w);
				assert(best != LLONG_MAX);
				tot += best;
				assert(nodes[a].dist == -inf);
			}
			assert(tot < 0);
		} else {
			assert(r.negCycle.empty());
			rep(i,0,n) {
				assert((bool)r.seen[i] == (nodes[i].dist != inf));
				if (r.seen[i]) assert(r.dist[i] == nodes[i].dist);
			}
			// parent chain reconstructs the distance
			rep(i,0,n) if (r.seen[i]) {
				int v = i; ll len = 0; int steps = 0;
				while (r.par[v] != -1) {
					int p = r.par[v]; ll best = LLONG_MAX;
					for (auto [x, w] : g[p]) if (x == v) best = min(best, w);
					len += best; v = p; assert(++steps <= n);
				}
				assert(v == s && len == r.dist[i]);
			}
		}
	}
	assert(withCycle > 1000);

	// difference constraints: x_v - x_u <= w, all vertices as sources
	rep(it,0,2000) {
		int n = rand() % 6 + 2, m = rand() % 10;
		vector<array<int, 3>> cons;
		vector<vector<pair<int, ll>>> g(n);
		rep(i,0,m) {
			int u = rand() % n, v = rand() % n, w = rand() % 11 - 5;
			cons.push_back({u, v, w});
			g[u].emplace_back(v, w);
		}
		vi all_(n); iota(all(all_), 0);
		auto r = spfa(g, all_);
		if (r.negCycle.empty()) {
			for (auto [u, v, w] : cons) assert(r.dist[v] - r.dist[u] <= w);
		} else {
			// brute force: no assignment in [-6,6]^n satisfies (n<=7 would be 13^7, too big) ->
			// instead verify the reported cycle is negative, which proves infeasibility
			ll tot = 0; auto& c = r.negCycle;
			rep(i,0,sz(c)) {
				int a = c[i], b = c[(i + 1) % sz(c)]; ll best = LLONG_MAX;
				for (auto [v, w] : g[a]) if (v == b) best = min(best, w);
				assert(best != LLONG_MAX); tot += best;
			}
			assert(tot < 0);
		}
	}
	// large: 2e5-node random graph with nonnegative weights (no cycle) and the
	// same graph plus one far-away negative cycle (must be found quickly)
	{
		int n = 200000;
		vector<vector<pair<int, ll>>> g(n);
		rep(i,0,n-1) g[i].emplace_back(i + 1, rand() % 100);
		rep(i,0,n) g[rand() % n].emplace_back(rand() % n, rand() % 100);
		auto r = spfa(g, {0});
		assert(r.negCycle.empty() && r.seen[n-1]);
		g[n-1].emplace_back(n - 2, -1000);
		g[n-2].emplace_back(n - 1, 5);
		auto r2 = spfa(g, {0});
		assert(sz(r2.negCycle) == 2);
	}
	cout << "Tests passed!" << endl;
}
