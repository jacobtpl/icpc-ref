#include "../utilities/template.h"

#include "../../content/graph/Dial.h"
#include "../../content/graph/Dijkstra.h"

void check(const vector<vector<pair<int, int>>>& g, vi src, const Dial& res) {
	int n = sz(g);
	auto& [dist, par] = res;
	auto r = dijkstra(g, src);
	rep(i,0,n) {
		assert((dist[i] != -1) == (bool)r.seen[i]);
		if (dist[i] == -1) continue;
		assert(dist[i] == r.dist[i]);
		// parent chain: ends at a source, consistent with edge weights
		int v = i, steps = 0;
		while (par[v] != -1) {
			int p = par[v];
			bool ok = false;
			for (auto [x, w] : g[p]) if (x == v && dist[p] + w == dist[v]) ok = true;
			assert(ok);
			v = p;
			assert(++steps <= n);
		}
		assert(dist[v] == 0 && find(all(src), v) != src.end());
	}
}

int main() {
	rep(it,0,5000) {
		int n = rand() % 9 + 1, m = rand() % 20, K = rand() % 5 + 1;
		if (rand() % 2) K = 1;
		vector<vector<pair<int, int>>> g(n);
		rep(i,0,m) {
			int a = rand() % n, b = rand() % n, w = rand() % (K + 1);
			g[a].emplace_back(b, w);
		}
		vi src;
		rep(i,0,rand() % 3 + 1) src.push_back(rand() % n);
		check(g, src, dial(g, src, K));
		if (K == 1) check(g, src, zeroOneBFS(g, src));
	}

	// grid: 0 = free move inside a row, 1 = row change; check a known answer
	{
		int R = 400, C = 500, n = R * C;
		vector<vector<pair<int, int>>> g(n);
		auto id = [&](int r, int c) { return r * C + c; };
		rep(r,0,R) rep(c,0,C) {
			if (c + 1 < C) g[id(r, c)].emplace_back(id(r, c + 1), 0),
			               g[id(r, c + 1)].emplace_back(id(r, c), 0);
			if (r + 1 < R) g[id(r, c)].emplace_back(id(r + 1, c), 1),
			               g[id(r + 1, c)].emplace_back(id(r, c), 1);
		}
		auto [dist, par] = zeroOneBFS(g, {id(0, 0)});
		rep(r,0,R) rep(c,0,C) assert(dist[id(r, c)] == r);
	}
	// larger K: 2e5-node random graph vs Dijkstra
	{
		int n = 200000, K = 10;
		vector<vector<pair<int, int>>> g(n);
		rep(i,0,n-1) g[i].emplace_back(i + 1, rand() % (K + 1));
		rep(i,0,n) g[rand() % n].emplace_back(rand() % n, rand() % (K + 1));
		auto res = dial(g, {0}, K);
		auto r = dijkstra(g, {0});
		rep(i,0,n) assert(res.dist[i] == r.dist[i]);
	}
	cout << "Tests passed!" << endl;
}
