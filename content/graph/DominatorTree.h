/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: Lengauer, Tarjan, "A Fast Algorithm for Finding
 * Dominators in a Flowgraph" (1979), simple version
 * Description: Dominator tree of a directed graph from root $r$.
 * $u$ dominates $v$ if every path from $r$ to $v$ passes through
 * $u$; the immediate dominator of $v \neq r$ is its closest
 * strict dominator, i.e. its parent in the dominator tree.
 * Returns idom for every vertex, with idom$[r] = r$ and $-1$ for
 * vertices unreachable from $r$. Vertices are 0-indexed;
 * self-loops and multi-edges are allowed. No recursion.
 * Usage: vi idom = domTree(g, r); // g[u] = out-neighbours of u
 * Time: O((V + E) \log V)
 * Status: stress-tested against brute force
 */
#pragma once

vi domTree(vector<vi>& g, int r) {
	int n = sz(g), k = 0;
	vi id(n, -1), ord, par(n), sd(n), lab(n), d(n), dom(n),
		res(n, -1), s;
	vector<vi> rg(n), bk(n);
	vector<pii> st{{r, -1}};
	while (!st.empty()) { // rg[v] = dfs numbers of preds of v
		auto [v, p] = st.back(); st.pop_back();
		if (p >= 0) rg[v].push_back(p);
		if (id[v] >= 0) continue;
		par[id[v] = k++] = p; ord.push_back(v);
		for (int u : g[v]) st.push_back({u, id[v]});
	}
	iota(all(sd), 0); lab = d = sd;
	auto eval = [&](int v) {
		for (s.clear(); d[v] != v; v = d[v]) s.push_back(v);
		for (int i = sz(s) - 1; i-- > 0;) {
			int x = s[i], y = s[i + 1];
			if (sd[lab[y]] < sd[lab[x]]) lab[x] = lab[y];
			d[x] = v;
		}
		return lab[s.empty() ? v : s[0]];
	};
	for (int w = k; --w > 0;) {
		for (int p : rg[ord[w]]) sd[w] = min(sd[w], sd[eval(p)]);
		bk[sd[w]].push_back(w);
		int p = d[w] = par[w];
		for (int v : bk[p]) {
			int u = eval(v);
			dom[v] = sd[u] < p ? u : p;
		}
		bk[p].clear();
	}
	res[r] = r;
	rep(w,1,k) {
		if (dom[w] != sd[w]) dom[w] = dom[dom[w]];
		res[ord[w]] = ord[dom[w]];
	}
	return res;
}
