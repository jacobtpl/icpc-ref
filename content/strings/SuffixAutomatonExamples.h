/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: https://cp-algorithms.com/string/suffix-automaton.html
 * Description: Common things to do with an SA after all strings
 * have been appended. Everything is 0-indexed.
 *  saOrder(sa): states sorted by increasing dis. This is a
 *   topological order of adj; iterate it backwards to process
 *   link-tree children before parents (e.g. to push counts up).
 *  saOcc(sa, s): cnt[v] = number of occurrences in s of each
 *   substring of state v (all equal). s must have been added.
 *   cnt[0] is meaningless.
 *  saDistinct(sa): number of distinct non-empty substrings over
 *   all added strings (fits in ll, can exceed int).
 *  saKth(sa, k): k-th (0 <= k < saDistinct) distinct non-empty
 *   substring in sorted order. Sorted by char, which is signed on
 *   most judges, so bytes >= 128 come first.
 *  saLcs(sa, t): {len, start} of a longest substring t[start,
 *   start+len) that is a substring of some added string.
 *   len = 0 if none.
 * Time: saOrder, saOcc, saDistinct O(n); saKth O(n + |ans| sigma);
 * saLcs O(|t| log sigma).
 * Usage:
 *  SA sa; int p = 0;
 *  for (char c : s) p = sa.append(p, c);
 *  vi cnt = saOcc(sa, s);          // cnt[p] == 1 (p = state of s)
 *  // # occurrences of w in s: walk w from state 0 as in
 *  // SuffixAutomaton.h to get v, answer is cnt[v] (0 if no path)
 *  ll D = saDistinct(sa);          // # distinct substrings of s
 *  string big = saKth(sa, D - 1);  // largest substring
 *  auto [len, st] = saLcs(sa, t);  // LCS of s, t is t.substr(st, len)
 *  // longest substring occurring >= 2 times in s:
 *  int best = 0;
 *  rep(v,1,sa.N) if (cnt[v] >= 2) best = max(best, sa.dis[v]);
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
