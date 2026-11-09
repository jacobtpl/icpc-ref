/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: https://cp-algorithms.com/string/suffix-automaton.html
 * Description: Common things to do with an \texttt{SA} once
 * all strings have been appended (all 0-indexed).
 * \texttt{saOrder}: states by increasing \texttt{dis} (a
 * topological order of \texttt{adj}; reversed, it visits
 * children before parents in the link tree).
 * \texttt{saOcc(s)}: for every state, the number of occurrences
 * in $s$ of each of its substrings ($s$ must have been added;
 * entry 0 is meaningless).
 * \texttt{saDistinct}: number of distinct non-empty substrings.
 * \texttt{saKth(k)}: $k$-th ($0 \le k <$ \texttt{saDistinct})
 * distinct non-empty substring in sorted order (by
 * \texttt{char}, which is usually signed).
 * \texttt{saLcs(t)}: \{length, start in $t$\} of a longest
 * substring of $t$ that is a substring of an added string.
 * Time: \texttt{saOrder}, \texttt{saOcc}, \texttt{saDistinct}
 * are $O(n)$, \texttt{saKth} is $O(n + |ans| \Sigma)$,
 * \texttt{saLcs} is $O(|t| \log \Sigma)$.
 * Status: stress-tested
 */
#pragma once

#include "SuffixAutomaton.h"

vi saOrder(SA& sa) {
	vi c(sa.N + 1), o(sa.N);
	rep(i,0,sa.N) c[sa.dis[i] + 1]++;
	rep(i,0,sa.N) c[i + 1] += c[i];
	rep(i,0,sa.N) o[c[sa.dis[i]]++] = i;
	return o;
}
vi saOcc(SA& sa, const string& s) {
	vi r(sa.N), o = saOrder(sa);
	int p = 0;
	for (char c : s) r[p = sa.adj[p][c]]++;
	for (int i = sa.N; --i;) r[sa.link[o[i]]] += r[o[i]];
	return r;
}
ll saDistinct(SA& sa) {
	ll r = 0;
	rep(i,1,sa.N) r += sa.dis[i] - sa.dis[sa.link[i]];
	return r;
}
string saKth(SA& sa, ll k) {
	vi o = saOrder(sa);
	vector<ll> d(sa.N, 1); // # of paths from each state, incl. empty
	for (int i = sa.N; i--;)
		for (auto [c, v] : sa.adj[o[i]]) d[o[i]] += d[v];
	string r;
	int p = 0;
	for (k++; k--;) for (auto [c, v] : sa.adj[p]) {
		if (k < d[v]) { r += c, p = v; break; }
		k -= d[v];
	}
	return r;
}
pii saLcs(SA& sa, const string& t) {
	pii r;
	int p = 0, l = 0;
	rep(i,0,sz(t)) {
		while (p && !sa.adj[p].count(t[i])) l = sa.dis[p = sa.link[p]];
		if (sa.adj[p].count(t[i])) p = sa.adj[p][t[i]], l++;
		r = max(r, pii(l, i + 1 - l));
	}
	return r;
}
