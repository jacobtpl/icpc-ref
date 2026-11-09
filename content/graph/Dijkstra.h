/**
 * Author: Devin
 * Date: 2026-11-08
 * License: CC0
 * Source: folklore
 * Description: Dijkstra over an adjacency list g[u] = {(v, w)} with any weight
 * type W. "What a distance is" lives in a small policy struct passed as the
 * last argument (default LLDist = plain ll). Result: dist[v] is valid iff
 * seen[v]; par[v] = predecessor on a shortest path (-1 for sources). Several
 * sources allowed, all starting at zero(). Requires add(d, w) to never be
 * less than d (nonnegative weights).
 *
 * Basic use:
 *   vector<vector<pair<int, ll>>> g(n);   // g[u].push_back({v, w});
 *   auto r = dijkstra(g, {s});
 *   if (r.seen[t]) cout << r.dist[t];     // path: t, r.par[t], r.par[r.par[t]], ..., s
 *
 * Custom distance: write a struct with exactly these 4 members and pass it.
 *   struct HopsDist {                       // lexicographic (length, #edges)
 *     typedef pair<ll, int> D;              // the distance type
 *     D zero() { return {0, 0}; }           // distance of a source
 *     D add(D d, ll w) { return {d.first + w, d.second + 1}; } // relax along w
 *     bool less(D a, D b) { return a < b; } // strict order, used by heap + relaxation
 *   };
 *   auto r = dijkstra(g, {s}, HopsDist{});  // r.dist[t].first, r.dist[t].second
 *
 * Huge numbers, e.g. CF 464E (w = 2^x, x <= 1e5): make D a handle into a
 * persistent structure, kept in the struct through a pointer. BigBits (a
 * persistent bitset with +2^x and compare) is in stress-tests/graph/DijkstraBigNum.cpp.
 *   struct Big {
 *     BigBits* t; typedef int D;                     // D = root of a persistent bitset
 *     D zero() { return 0; }
 *     D add(D d, int x) { return t->add(d, x); }     // d + 2^x, new root
 *     bool less(D a, D b) { return t->less(a, b); }  // compare by hash descent
 *   };
 *   BigBits bb; auto r = dijkstra(g, {s}, Big{&bb}); // g: vector<vector<pair<int,int>>>
 *
 * Implicit graph (edges generated on the fly): edges(u, f) must call f(v, w)
 * for every edge u -> v.
 *   auto r = dijkstraGen(n, {s}, LLDist{}, [&](int u, auto f) {
 *     if (u + 1 < n) f(u + 1, 1);
 *     if (2 * u < n) f(2 * u, 3);
 *   });
 * Time: O((V+E) log E)
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
