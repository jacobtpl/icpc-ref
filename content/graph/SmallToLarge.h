/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: Folklore (DSU on tree / sack), own non-recursive formulation
 * Description: Small-to-large on a rooted tree. For every vertex $v$
 * reachable from $r$, calls ans($v$) at a moment when exactly the
 * vertices of the subtree of $v$ are added: the data of the heavy
 * child is kept, light subtrees are re-added. g is a 0-indexed
 * adjacency list of the tree, directed away from $r$ or undirected,
 * $N \ge 1$. No recursion. Each vertex is added (and removed) at
 * most $1 + \log_2 N$ times; hooks are called in no particular
 * order, so they must commute. Everything is removed at the end.
 * Usage:
 *  vi cnt(C), res(n); int d = 0; // distinct colours in subtree
 *  sack(g, 0, [\&](int v) { d += !cnt[col[v]]++; }, [\&](int v)
 *   { d -= !--cnt[col[v]]; }, [\&](int v) { res[v] = d; });
 * Time: $O(N \log N)$ hook calls + $O(N)$
 * Status: stress-tested
 */
#pragma once

template<class A, class R, class F>
void sack(vector<vi>& g, int r, A add, R rem, F ans) {
	int n = sz(g);
	vi p(n, -1), s(n, 1), q{r}, in(n), at(n);
	rep(i,0,sz(q)) for (int y : g[q[i]]) if (y != p[q[i]])
		p[y] = q[i], q.push_back(y);
	for (int i = sz(q); --i;) s[p[q[i]]] += s[q[i]];
	for (int v : q) { // preorder, heavy child first
		int h = -1, t = in[v] + 1;
		at[t - 1] = v;
		for (int y : g[v]) if (y != p[v] && (h < 0 || s[y] > s[h]))
			h = y;
		if (h >= 0) in[h] = t, t += s[h];
		for (int y : g[v]) if (y != p[v] && y != h)
			in[y] = t, t += s[y];
	}
	for (int i = sz(q); i--;) {
		int v = at[i];
		add(v);
		if (s[v] > 1) rep(j,i+1+s[at[i+1]],i+s[v]) add(at[j]);
		ans(v);
		if (!i || p[v] != at[i-1]) rep(j,i,i+s[v]) rem(at[j]);
	}
}
