/**
 * Author: Devin
 * Date: 2026-11-21
 * License: CC0
 * Description: Max flow (Dinic). addEdge returns an edge id; after calc,
 *  flow(id) is the flow on that edge (negative = backwards, if rcap > 0).
 *  For an undirected edge pass rcap = c. calc pushes at most lim more flow
 *  and can be called again, also after adding edges.
 *  leftOfMinCut is valid after a calc that was not stopped by lim.
 *  Change F for another capacity type (not floating point).
 * Usage: Dinic d(n); int id = d.addEdge(a, b, cap); ll f = d.calc(s, t);
 * Time: O(V^2 E); O(E sqrt(V)) for unit capacities / bipartite matching.
 *  dfs recurses once per path vertex: run with ulimit -s unlimited locally.
 * Status: stress-tested against PushRelabel and EdmondsKarp
 */
#pragma once

struct Dinic {
	using F = ll;
	struct Edge { int to, rev; F c, cap; }; // c = residual capacity
	vector<vector<Edge>> adj;
	vector<pii> pos; // edge id -> (a, index in adj[a])
	vi lvl, ptr, q;
	Dinic(int n) : adj(n), lvl(n), ptr(n), q(n) {}
	int addEdge(int a, int b, F c, F rcap = 0) {
		pos.push_back({a, sz(adj[a])});
		adj[a].push_back({b, sz(adj[b]) + (a == b), c, c});
		adj[b].push_back({a, sz(adj[a]) - 1, rcap, rcap});
		return sz(pos) - 1;
	}
	F flow(int id) {
		Edge& e = adj[pos[id].first][pos[id].second];
		return e.cap - e.c;
	}
	F dfs(int v, int s, F f) { // walks back from t, so only vertices reachable from s are visited
		if (v == s) return f;
		F left = f;
		for (int& i = ptr[v]; i < sz(adj[v]); i++) {
			Edge& e = adj[v][i], &r = adj[e.to][e.rev];
			if (r.c && lvl[e.to] && lvl[e.to] == lvl[v] - 1) {
				F p = dfs(e.to, s, min(left, r.c));
				r.c -= p, e.c += p;
				if (!(left -= p)) return f;
			}
		}
		return f - left;
	}
	F calc(int s, int t, F lim = numeric_limits<F>::max()) {
		assert(s != t);
		F flow = 0;
		while (flow < lim) {
			fill(all(lvl), 0), fill(all(ptr), 0);
			int qi = 0, qe = 0;
			q[qe++] = s, lvl[s] = 1;
			while (qi < qe && !lvl[t]) {
				int v = q[qi++];
				for (Edge& e : adj[v]) if (e.c && !lvl[e.to])
					q[qe++] = e.to, lvl[e.to] = lvl[v] + 1;
			}
			if (!lvl[t]) break;
			flow += dfs(t, s, lim - flow);
		}
		return flow;
	}
	bool leftOfMinCut(int v) { return lvl[v] != 0; }
};
