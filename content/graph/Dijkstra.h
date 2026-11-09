/**
 * Author: Devin
 * Date: 2026-11-08
 * License: CC0
 * Source: folklore
 * Description: Dijkstra with a pluggable distance type $D$ and comparator.
 * \texttt{edges(u, du, push)} is called exactly once per node $u$, when $u$ is
 * finalized with distance \texttt{du}; it must call \texttt{push(v, dv)} for every
 * edge $u \to v$ with the candidate distance \texttt{dv}. $D$ only needs the strict
 * weak order \texttt{less}; it can be \texttt{ll}, a pair, a float, or a handle into
 * a persistent structure representing a huge number (e.g. CF 464E, where edge
 * weights are $2^x$). Edges may be generated on the fly (implicit graphs).
 * Several sources with arbitrary initial distances are allowed.
 * \texttt{dist[v]} is valid iff \texttt{seen[v]}; \texttt{par[v]} is the predecessor
 * on a shortest path (-1 for sources). \texttt{dijkstra(g, s)} is the plain
 * adjacency-list version (\texttt{g[u] = \{(v, w)\}}, unreachable = inf).
 * Usage: auto r = dijkstra<ll>(n, {{s, 0}}, less<ll>(),
 *   [\&](int u, ll du, auto push) { for (auto [v, w] : g[u]) push(v, du + w); });
 * Time: O((V+E) \log E)
 * Status: stress-tested against Bellman-Ford; used to solve CF 464E
 */
#pragma once

template<class D> struct DijkstraRes {
	vector<D> dist; vi par; vector<char> seen;
};
template<class D, class Less, class Edges>
DijkstraRes<D> dijkstra(int n, vector<pair<int, D>> src,
                        Less less, Edges edges) {
	DijkstraRes<D> r{vector<D>(n), vi(n, -1), vector<char>(n)};
	vector<char> done(n);
	auto cmp = [&](const pair<D, int>& a, const pair<D, int>& b) {
		return less(b.first, a.first);
	};
	priority_queue<pair<D, int>, vector<pair<D, int>>,
	               decltype(cmp)> pq(cmp);
	auto push = [&](int v, D dv, int p) {
		if (r.seen[v] && !less(dv, r.dist[v])) return;
		r.dist[v] = dv; r.par[v] = p; r.seen[v] = 1;
		pq.emplace(dv, v);
	};
	for (auto& [s, d0] : src) push(s, d0, -1);
	while (!pq.empty()) {
		auto [du, u] = pq.top(); pq.pop();
		if (done[u]) continue;
		done[u] = 1;
		edges(u, du, [&](int v, D dv) { push(v, dv, u); });
	}
	return r;
}

vector<ll> dijkstra(const vector<vector<pair<int, ll>>>& g, int s) {
	auto r = dijkstra<ll>(sz(g), {{s, 0}}, less<ll>(),
		[&](int u, ll du, auto push) {
			for (auto [v, w] : g[u]) push(v, du + w);
		});
	rep(i,0,sz(g)) if (!r.seen[i]) r.dist[i] = LLONG_MAX;
	return r.dist;
}
