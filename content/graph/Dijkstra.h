/**
 * Author: Devin
 * Date: 2026-11-08
 * License: CC0
 * Source: folklore
 * Description: Dijkstra over \texttt{g[u] = \{(v, w)\}} with any weight type $W$.
 * What a "distance" is lives in a small policy struct passed as the last
 * argument; to use a custom distance, write your own with these 4 members:
 *   struct MyDist {
 *     typedef ... D;         // distance type
 *     D zero();              // distance of a source
 *     D add(D d, W w);       // distance after taking an edge of weight w
 *     bool less(D a, D b);   // strict order: used for both heap and relaxation
 *   };
 * Default \texttt{LLDist} is plain \texttt{ll}. Keep state (e.g. a persistent
 * tree) in the struct as a pointer. Example, CF 464E (w = $2^x$, $x \le 10^5$;
 * distances are roots of a persistent bitset, see DijkstraBigNum stress test):
 *   struct Big { BigBits* t; typedef int D;
 *     D zero() { return 0; }
 *     D add(D d, int x) { return t->add(d, x); }
 *     bool less(D a, D b) { return t->less(a, b); } };
 *   auto r = dijkstra(g, {s}, Big{\&bb});
 * Other overrides: \texttt{pair<ll,int>} for (length, hops), double, etc.
 * All sources start at \texttt{zero()}. \texttt{dist[v]} is valid iff
 * \texttt{seen[v]}; \texttt{par[v]} = predecessor (-1 for sources).
 * \texttt{dijkstraGen(n, src, dist, edges)} is the same on an implicit graph:
 * \texttt{edges(u, f)} must call \texttt{f(v, w)} for every edge $u \to v$.
 * Usage: auto r = dijkstra(g, {s}); // r.dist[v] is ll
 * Time: O((V+E) \log E)
 * Status: stress-tested against Bellman-Ford; solves CF 464E
 */
#pragma once

struct LLDist {
	typedef ll D;
	D zero() { return 0; }
	D add(D d, ll w) { return d + w; }
	bool less(D a, D b) { return a < b; }
};
template<class D> struct DijkstraRes {
	vector<D> dist; vi par; vector<char> seen;
};
template<class Dist, class Edges>
DijkstraRes<typename Dist::D> dijkstraGen(int n, vi src, Dist dist,
                                          Edges edges) {
	typedef typename Dist::D D;
	DijkstraRes<D> r{vector<D>(n), vi(n, -1), vector<char>(n)};
	vector<char> done(n);
	auto cmp = [&](const pair<D, int>& a, const pair<D, int>& b) {
		return dist.less(b.first, a.first);
	};
	priority_queue<pair<D, int>, vector<pair<D, int>>,
	               decltype(cmp)> pq(cmp);
	auto relax = [&](int v, D dv, int p) {
		if (r.seen[v] && !dist.less(dv, r.dist[v])) return;
		r.dist[v] = dv; r.par[v] = p; r.seen[v] = 1;
		pq.emplace(dv, v);
	};
	for (int s : src) relax(s, dist.zero(), -1);
	while (!pq.empty()) {
		auto [du, u] = pq.top(); pq.pop();
		if (done[u]) continue;
		done[u] = 1;
		edges(u, [&](int v, auto w) { relax(v, dist.add(du, w), u); });
	}
	return r;
}
template<class W, class Dist = LLDist>
DijkstraRes<typename Dist::D> dijkstra(
		const vector<vector<pair<int, W>>>& g, vi src, Dist dist = {}) {
	return dijkstraGen(sz(g), src, dist, [&](int u, auto f) {
		for (auto& [v, w] : g[u]) f(v, w);
	});
}
