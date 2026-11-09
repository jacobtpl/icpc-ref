#include "../utilities/template.h"

#include "../../content/graph/ZeroOneBFS.h"
#include "../../content/graph/Dijkstra.h"

int main() {
	rep(it,0,5000) {
		int n = rand() % 9 + 1, m = rand() % 20;
		vector<vector<pair<int, int>>> g(n);
		vector<vector<pair<int, ll>>> gl(n);
		rep(i,0,m) {
			int a = rand() % n, b = rand() % n, w = rand() % 2;
			g[a].emplace_back(b, w);
			gl[a].emplace_back(b, w);
		}
		int k = rand() % 3 + 1;
		vi src;
		rep(i,0,k) src.push_back(rand() % n);

		auto [dist, par] = zeroOneBFS(g, src);
		auto r = dijkstra<ll>(n, [&] {
			vector<pair<int, ll>> s;
			for (int x : src) s.emplace_back(x, 0);
			return s;
		}(), less<ll>(), [&](int u, ll du, auto push) {
			for (auto [v, w] : gl[u]) push(v, du + w);
		});
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
	cout << "Tests passed!" << endl;
}
