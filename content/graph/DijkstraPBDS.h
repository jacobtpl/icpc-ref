/**
 * Author: Devin
 * Date: 2026-11-19
 * License: CC0
 * Source: folklore
 * Description: Same as Dijkstra.h (same arguments, policy structs and
 * DijkstraRes), but the heap is a pbds pairing heap holding each vertex at most
 * once: an improved distance does modify (decrease-key) instead of pushing a
 * second copy, so the heap has at most V entries instead of up to E. Use it
 * when Dist::less is expensive (big numbers, CF 464E) or many relaxations
 * succeed: about 3 comparisons per edge there vs about 40 for Dijkstra.h.
 * On random graphs with ll distances Dijkstra.h is up to 2x faster.
 *   auto r = dijkstraPBDS(g, {s});                // same as dijkstra(g, {s})
 *   auto r = dijkstraPBDSGen(n, {s}, Big{&bb}, edges);
 * thin_heap_tag also works as the last template argument; rc_binomial_heap_tag
 * gives wrong answers with modify.
 * Time: O(E log V)
 * Status: stress-tested against Dijkstra.h and Bellman-Ford
 */
#pragma once

#include <ext/pb_ds/priority_queue.hpp>
#include "Dijkstra.h"

template<class Dist, class Edges, class Tag = __gnu_pbds::pairing_heap_tag>
DijkstraRes<typename Dist::D> dijkstraPBDSGen(int n, vi src, Dist dist,
                                              Edges edges) {
	typedef typename Dist::D D; typedef pair<D, int> P;
	DijkstraRes<D> r{vector<D>(n), vi(n, -1), vector<char>(n)};
	auto cmp = [&](const P& a, const P& b) {
		return dist.less(b.first, a.first);
	};
	typedef __gnu_pbds::priority_queue<P, decltype(cmp), Tag> PQ;
	PQ pq(cmp);
	vector<typename PQ::point_iterator> it(n);
	auto relax = [&](int v, D dv, int p) {
		if (r.seen[v] && !dist.less(dv, r.dist[v])) return;
		r.dist[v] = dv; r.par[v] = p;
		if (r.seen[v]) pq.modify(it[v], {dv, v});
		else r.seen[v] = 1, it[v] = pq.push({dv, v});
	};
	for (int s : src) relax(s, dist.zero(), -1);
	while (!pq.empty()) {
		auto [du, u] = pq.top(); pq.pop();
		edges(u, [&](int v, auto w) { relax(v, dist.add(du, w), u); });
	}
	return r;
}
template<class W, class Dist = LLDist>
DijkstraRes<typename Dist::D> dijkstraPBDS(
		const vector<vector<pair<int, W>>>& g, vi src, Dist dist = {}) {
	return dijkstraPBDSGen(sz(g), src, dist, [&](int u, auto f) {
		for (auto& [v, w] : g[u]) f(v, w);
	});
}
