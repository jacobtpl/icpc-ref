/**
 * Author: Devin
 * Date: 2026-11-08
 * License: CC0
 * Source: folklore
 * Description: SPFA (queue-based Bellman-Ford) on \texttt{g[u] = \{(v, w)\}}
 * with negative weights allowed; several sources (all at distance 0).
 * Returns \texttt{dist/par} (valid iff \texttt{seen[v]}; unreachable otherwise)
 * and \texttt{negCycle}: the vertices of a negative cycle reachable from a
 * source, in order, or empty if none (then all dists are final).
 * When to use: negative edges but no cycles expected (difference constraints,
 * potentials for min-cost flow / Johnson). Difference constraints
 * $x_v - x_u \le w$: add edge $u \to v$ weight $w$, run from all vertices as
 * sources (or a super-source), feasible iff \texttt{negCycle} empty, then
 * $x_v = \texttt{dist[v]}$ is a solution ($x_v \ge x_u + w$ becomes $v \to u$
 * with weight $-w$). Worst case is $O(VE)$ and it is easy to hit: a long
 * path plus random shortcut edges already takes minutes at $V = 2 \cdot 10^5$,
 * and CF hacks target it. For nonnegative weights use Dijkstra; for a
 * hack-proof negative-weight run use \texttt{BellmanFord.h}.
 * Usage: auto r = spfa(g, {s}); if (!r.negCycle.empty()) ...; else r.dist[v]
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
