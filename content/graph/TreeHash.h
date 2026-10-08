/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: own work; Aho, Hopcroft, Ullman (1974) tree isomorphism
 * Description: Deterministic tree isomorphism. treeId roots the
 * tree containing $r$ at $r$ and returns, for every vertex, an id
 * of its rooted subtree: two ids are equal iff the subtrees are
 * isomorphic as rooted trees, also across calls and across
 * different trees (ids live in the global map tid). unrootedId
 * returns the smaller root id over the (one or two) centroids; two
 * trees are isomorphic iff these are equal. Only compare ids of
 * the same kind. g is a 0-indexed symmetric adjacency list of a
 * tree with $N \ge 1$ vertices. No recursion.
 * Time: O(N \log N)
 * Usage: vi id = treeId(g, 0); // id[v] == id[u]: same subtree
 *  bool iso = unrootedId(g1) == unrootedId(g2);
 * Status: stress-tested
 */
#pragma once

map<vi, int> tid;
vi siz; // subtree sizes from the last treeId call
vi treeId(vector<vi>& g, int r) {
	int n = sz(g);
	vi id(n), par(n, -1), q{r};
	siz.assign(n, 1);
	rep(i,0,sz(q)) for (int u : g[q[i]]) if (u != par[q[i]])
		par[u] = q[i], q.push_back(u);
	for (int i = sz(q); i--;) {
		int v = q[i]; vi c;
		for (int u : g[v]) if (u != par[v])
			c.push_back(id[u]), siz[v] += siz[u];
		sort(all(c));
		id[v] = tid.emplace(c, sz(tid)).first->second;
	}
	return id;
}
int unrootedId(vector<vi>& g) {
	int n = sz(g), res = INT_MAX;
	treeId(g, 0);
	vi s = siz;
	rep(v,0,n) {
		int m = n - s[v];
		for (int u : g[v]) if (s[u] < s[v]) m = max(m, s[u]);
		if (2 * m <= n) res = min(res, treeId(g, v)[v]);
	}
	return res;
}
