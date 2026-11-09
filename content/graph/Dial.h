/**
 * Author: Devin
 * Date: 2026-11-08
 * License: CC0
 * Source: folklore
 * Description: Dial's algorithm: single/multi-source shortest paths when all
 * edge weights are integers in $[0, K]$ (\texttt{g[u] = \{(v, w)\}}). Bucket
 * queue of $K+1$ vectors indexed by $d \bmod (K+1)$ replaces the heap.
 * \texttt{zeroOneBFS} is the $K=1$ case (weight-0 edges land in the current
 * bucket and are processed before moving on); it is as fast as the usual
 * deque 0-1 BFS, so there is only one implementation. Unreachable = -1;
 * \texttt{par[v]} is the predecessor (-1 for sources).
 * Usage: auto [dist, par] = zeroOneBFS(g, {s});
 *        auto [dist, par] = dial(g, {s1, s2}, K);
 * Time: O(VK + E)
 * Status: stress-tested against Dijkstra
 */
#pragma once

struct Dial { vi dist, par; };
Dial dial(const vector<vector<pair<int, int>>>& g, vi src, int K) {
	int n = sz(g), left = 0;
	vi dist(n, -1), par(n, -1);
	vector<vi> b(K + 1);
	for (int s : src) if (dist[s] == -1)
		dist[s] = 0, b[0].push_back(s), left++;
	for (int d = 0; left; d++) {
		vi& cur = b[d % (K + 1)];
		while (!cur.empty()) {
			int u = cur.back(); cur.pop_back(); left--;
			if (dist[u] != d) continue;
			for (auto [v, w] : g[u]) {
				int nd = d + w;
				if (dist[v] == -1 || nd < dist[v]) {
					dist[v] = nd; par[v] = u;
					b[nd % (K + 1)].push_back(v); left++;
				}
			}
		}
	}
	return {dist, par};
}
Dial zeroOneBFS(const vector<vector<pair<int, int>>>& g, vi src) {
	return dial(g, src, 1);
}
