/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: own work; folklore prefix/suffix rerooting technique
 * Description: Rerooting DP on a forest given as 0-indexed symmetric
 *  adjacency lists $g$ (no self-loops or multi-edges). Returns, for every
 *  vertex $r$, $f(x, r, -1)$ where $x$ is the merge, in the order of
 *  $g[r]$, of the values sent by the neighbours of $r$. The value $v$
 *  sends to $g[v][i]$ is $f(x, v, i)$, $x$ being the merge of the values
 *  sent to $v$ by its other neighbours. $mg$ must be associative with
 *  identity $e$; it needs no inverse and need not be commutative.
 *  Non-recursive.
 * Usage:
 *  // longest path from each vertex, w[v][i] = weight of g[v][i]
 *  auto mx = [](ll a, ll b) { return max(a, b); };
 *  vector<ll> far = reroot(g, 0LL, mx,
 *   [\&](ll x, int v, int i) { return i < 0 ? x : x + w[v][i]; });
 * Time: O(N) calls to $mg$ and $f$
 * Status: stress-tested
 */
#pragma once

template<class T, class M, class F>
vector<T> reroot(const vector<vi>& g, T e, M mg, F f) {
	int n = sz(g);
	vi p(n, -1), q;
	rep(r,0,n) if (p[r] < 0) {
		p[r] = r; q.pb(r);
		for (int j = sz(q) - 1; j < sz(q); j++)
			for (int u : g[q[j]]) if (p[u] < 0)
				p[u] = q[j], q.pb(u);
	}
	vector<T> lo(n, e), hi(n, e), res(n, e), pre;
	for (int j = n; j--;) {
		int v = q[j], k = -1;
		T x = e;
		rep(i,0,sz(g[v]))
			if (g[v][i] == p[v]) k = i;
			else x = mg(x, lo[g[v][i]]);
		if (k >= 0) lo[v] = f(x, v, k);
	}
	for (int v : q) {
		int d = sz(g[v]);
		pre.assign(d + 1, e);
		rep(i,0,d) {
			int u = g[v][i];
			pre[i+1] = mg(pre[i], u == p[v] ? hi[v] : lo[u]);
		}
		res[v] = f(pre[d], v, -1);
		T s = e;
		for (int i = d; i--;) {
			int u = g[v][i];
			if (u != p[v]) hi[u] = f(mg(pre[i], s), v, i);
			s = mg(u == p[v] ? hi[v] : lo[u], s);
		}
	}
	return res;
}
