/**
 * Author: Devin
 * Date: 2026-11-08
 * License: CC0
 * Source: folklore
 * Description: Dial's algorithm: shortest paths when all edge weights are
 * integers in [0, K]; g[u] = {(v, w)}. A bucket queue of K+1 vectors indexed
 * by d % (K+1) replaces the heap. zeroOneBFS is the K = 1 case (weight-0
 * edges go into the current bucket, so they are processed before moving on)
 * and is as fast as the usual deque 0-1 BFS, hence a single implementation.
 * Several sources allowed. Unreachable = -1; par[v] = predecessor (-1 for sources).
 *
 * Example: grid, moving along a row is free, changing row costs 1.
 *   vector<vector<pair<int, int>>> g(R * C);
 *   auto id = [&](int r, int c) { return r * C + c; };
 *   g[id(r, c)].push_back({id(r, c + 1), 0});   // ... for all cells / directions
 *   g[id(r, c)].push_back({id(r + 1, c), 1});
 *   auto [dist, par] = zeroOneBFS(g, {id(0, 0)});
 *   for (int v = t; v != -1; v = par[v]) ...;   // path from t back to the source
 *
 *   auto [dist, par] = dial(g, {s1, s2}, K);   // weights in [0, K], multi-source
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
