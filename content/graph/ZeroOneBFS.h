/**
 * Author: Devin
 * Date: 2026-11-08
 * License: CC0
 * Source: folklore
 * Description: 0-1 BFS: single/multi-source shortest paths when every edge
 * weight is 0 or 1 (\texttt{g[u] = \{(v, w)\}}, $w \in \{0,1\}$). Deque replaces
 * the heap: weight-0 edges push to the front, weight-1 to the back.
 * Unreachable = -1; \texttt{par[v]} is the predecessor (-1 for sources).
 * Typical uses: grid "number of turns / walls broken" problems, graph with
 * free edges inside components. For small integer weights $0..K$ use a bucket
 * queue (Dial's algorithm, $O(VK+E)$) or just Dijkstra.
 * Usage: auto [dist, par] = zeroOneBFS(g, {s});
 * Time: O(V + E)
 * Status: stress-tested against Dijkstra
 */
#pragma once

struct ZeroOneBFS { vi dist, par; };
ZeroOneBFS zeroOneBFS(const vector<vector<pair<int, int>>>& g, vi src) {
	int n = sz(g);
	vi dist(n, -1), par(n, -1);
	deque<int> dq;
	for (int s : src) if (dist[s] == -1) dist[s] = 0, dq.push_back(s);
	while (!dq.empty()) {
		int u = dq.front(); dq.pop_front();
		for (auto [v, w] : g[u]) {
			int nd = dist[u] + w;
			if (dist[v] == -1 || nd < dist[v]) {
				dist[v] = nd; par[v] = u;
				if (w) dq.push_back(v); else dq.push_front(v);
			}
		}
	}
	return {dist, par};
}
