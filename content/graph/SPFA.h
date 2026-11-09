/**
 * Author: Devin
 * Date: 2026-11-08
 * License: CC0
 * Source: folklore
 * Description: SPFA (queue-based Bellman-Ford) on g[u] = {(v, w)}; negative
 * weights allowed, several sources (all at distance 0). Result: dist/par are
 * valid iff seen[v] (par = -1 for sources). negCycle = vertices of a negative
 * cycle reachable from a source, in edge order (c0 -> c1 -> ... -> ck -> c0);
 * empty means no such cycle and all dists are final.
 * When to use: negative edges but (usually) no negative cycle: difference
 * constraints, potentials for Johnson / min-cost flow. Worst case is O(VE)
 * and easy to hit (a long path plus random shortcut edges already takes
 * minutes at V = 2e5; CF hacks target it). Nonnegative weights -> Dijkstra;
 * hack-proof negative weights -> BellmanFord.h.
 * Detection: a relaxation in pass n-1 proves a negative cycle (Bellman-Ford
 * bound); the cycle is read off the parent graph, which is also checked every
 * n relaxations so an early cycle is found without waiting n passes.
 *
 * Basic use:
 *   vector<vector<pair<int, ll>>> g(n);   // g[u].push_back({v, w});
 *   auto r = spfa(g, {s});
 *   if (!r.negCycle.empty()) { ... }      // e.g. print r.negCycle
 *   else if (r.seen[t]) cout << r.dist[t];
 *
 * Difference constraints x[v] - x[u] <= w. (x[v] - x[u] >= w is
 * x[u] - x[v] <= -w; equality = both.) Every vertex is a source, which is
 * the same as a super-source with 0-edges to all:
 *   for (auto [u, v, w] : cons) g[u].push_back({v, w});
 *   vi allv(n); iota(all(allv), 0);
 *   auto r = spfa(g, allv);
 *   if (!r.negCycle.empty()) { ... }      // infeasible
 *   else { ... }                          // x[v] = r.dist[v] (all <= 0) is a solution
 * Time: O(VE) worst case, usually much faster
 * Status: stress-tested against BellmanFord
 */
#pragma once

struct SPFA { vector<ll> dist; vi par; vector<char> seen; vi negCycle; };
SPFA spfa(const vector<vector<pair<int, ll>>>& g, vi src) {
	int n = sz(g), pass = 0, cnt = 0, stamp = 0;
	SPFA r{vector<ll>(n), vi(n, -1), vector<char>(n), {}};
	vi vis(n); vector<char> inq(n);
	queue<int> q;
	for (int s : src) if (!r.seen[s])
		r.seen[s] = inq[s] = 1, q.push(s);
	q.push(-1); // pass separator
	while (sz(q) > 1) {
		int u = q.front(); q.pop();
		if (u == -1) { pass++; q.push(-1); continue; }
		inq[u] = 0;
		for (auto [v, w] : g[u]) {
			ll d = r.dist[u] + w;
			if (r.seen[v] && d >= r.dist[v]) continue;
			r.dist[v] = d; r.par[v] = u; r.seen[v] = 1;
			// a relaxation in pass n-1 proves a negative cycle; additionally
			// look for one in the parent graph every n relaxations (cheap).
			if (pass >= n - 1 || ++cnt % n == 0) {
				int x = v; stamp++;
				for (; x != -1 && vis[x] != stamp; x = r.par[x]) vis[x] = stamp;
				if (x != -1) {
					int y = x;
					do { r.negCycle.push_back(y); y = r.par[y]; } while (y != x);
					reverse(all(r.negCycle));
					return r;
				}
			}
			if (!inq[v]) inq[v] = 1, q.push(v);
		}
	}
	return r;
}
