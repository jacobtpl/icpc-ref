/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: Schrijver, Combinatorial Optimization, ch. 41;
 *  https://codeforces.com/blog/entry/69287
 * Description: Largest set of elements $0..n-1$ independent in two
 *  matroids, by shortest augmenting paths in the exchange graph.
 *  A matroid is any struct with \texttt{init(S)}, which loads an
 *  independent set $S$ (list of elements), and \texttt{ok(e)}, which
 *  tells if $S \cup \{e\}$ is independent for $e \notin S$.
 *  \texttt{Graphic}: element $i$ is edge \texttt{ed[i]} of a graph
 *  on vertices $0..V-1$, independent = forest (loops, multi-edges
 *  fine). \texttt{Colorful}: element $i$ has colour \texttt{col[i]}
 *  in $[0,$ \texttt{sz(cap)}$)$, at most \texttt{cap[c]} $\ge 0$ of
 *  colour $c$. Starts from a greedy common independent set.
 *  Returns the elements of the answer in increasing order.
 *  The exchange graph takes $O(rn)$ memory.
 * Time: With answer size $r$: $O(r^2)$ calls to init and
 *  $O(r^2 n)$ calls to ok. For the matroids below
 *  $O(r^2 (n + V))$ times $\alpha$. Rounds = $r$ minus greedy size:
 *  random $n=5000$, $r=500$ in 0.3s, $r=1000$ in 4s.
 * Usage: Graphic a{V, edges}; Colorful b{col, cap};
 *  vi res = matroidIsect(sz(edges), a, b);
 * Status: stress-tested against brute force
 */
#pragma once

#include "../data-structures/UnionFind.h"

struct Graphic {
	int V; vector<pii> ed; UF uf{0};
	void init(const vi& S) {
		uf = UF(V);
		for (int i : S) uf.join(ed[i].first, ed[i].second);
	}
	bool ok(int i) {
		return !uf.sameSet(ed[i].first, ed[i].second);
	}
};
struct Colorful {
	vi col, cap, c;
	void init(const vi& S) {
		c = cap;
		for (int i : S) c[col[i]]--;
	}
	bool ok(int i) { return c[col[i]] > 0; }
};

template<class A, class B>
vi matroidIsect(int n, A& a, B& b) {
	vector<bool> in(n);
	vi S; a.init(S); b.init(S);
	rep(i,0,n) if (a.ok(i) && b.ok(i)) // greedy start
		in[i] = 1, S.pb(i), a.init(S), b.init(S);
	for (;;) {
		vi I, J, q, par(n, -2), snk(n);
		vector<vi> g(n);
		rep(i,0,n) if (in[i]) I.pb(i);
		a.init(I); b.init(I);
		rep(i,0,n) if (!in[i]) {
			if (a.ok(i)) par[i] = -1, q.pb(i);
			snk[i] = b.ok(i);
		}
		for (int y : I) { // edges y->x (a), x->y (b)
			J.clear();
			for (int z : I) if (z != y) J.pb(z);
			a.init(J); b.init(J);
			rep(x,0,n) if (!in[x]) {
				if (a.ok(x)) g[y].pb(x);
				if (b.ok(x)) g[x].pb(y);
			}
		}
		int t = -1;
		rep(i,0,sz(q)) { // BFS: path must be shortest
			int v = q[i];
			if (snk[v]) { t = v; break; }
			for (int w : g[v]) if (par[w] == -2)
				par[w] = v, q.pb(w);
		}
		if (t < 0) return I;
		for (; t >= 0; t = par[t]) in[t] = !in[t];
	}
}
