/**
 * Author: Devin
 * Date: 2026-11-17
 * License: CC0
 * Source: folklore; UCup 5 Stage 1 C
 * Description: Boruvka MST skeleton for graphs too big to list (complete
 *  graphs with weights given by a formula or a data structure).
 *  Each round calls f(comp, best): comp[v] is the component id of v (a
 *  vertex index); f must set best[c] = {w, u, v} to a minimum-weight edge
 *  with comp[u] == c != comp[v], for every component id c that has one
 *  (best is pre-filled with v = -1). At most log2(N) rounds.
 *  Returns {total weight, edges}; the result is a forest with
 *  sz(edges) < N-1 if the graph is disconnected.
 *  Typical f: sweep / segment tree keeping the best and the best from a
 *  different component (top two by component) for every query.
 * Time: O(log N) calls of f
 * Status: stress-tested against Kruskal
 */
#pragma once

#include "../data-structures/UnionFind.h"

template<class W, class F>
pair<W, vector<pii>> boruvka(int n, F f) {
	UF uf(n); W tot{}; vector<pii> ed; vi comp(n);
	for (bool any = 1; any;) {
		any = 0;
		rep(i,0,n) comp[i] = uf.find(i);
		vector<tuple<W, int, int>> best(n, {W{}, -1, -1});
		f(comp, best);
		for (auto [w, u, v] : best) if (v >= 0 && uf.join(u, v))
			tot += w, ed.push_back({u, v}), any = 1;
	}
	return {tot, ed};
}
