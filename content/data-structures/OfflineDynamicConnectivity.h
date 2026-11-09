/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: folklore (segment tree over time, "offline dynamic connectivity");
 * written from scratch on top of UnionFindRollback.h
 * Description: Offline dynamic connectivity. Call add/rem/query in
 * chronological order, then solve(f) calls f(i) for each query
 * $i = 0, 1, \dots$ (the value returned by query()) in order, with uf
 * holding exactly the edges present at that moment.
 * Nodes are 0-indexed. Multi-edges and self-loops are fine: rem(a, b)
 * removes one copy of the undirected edge and requires that one is present.
 * Edges never removed stay until the end. Call solve at most once.
 * The number of components is $N - sz(uf.st)/2$.
 * Time: $O(E \log Q \log N + Q \log N)$ for $E$ adds and $Q$ queries.
 * Memory: $O(E \log Q)$
 * Usage: DynCon d(n); d.add(a, b); qs.push_back({u, v}); d.query();
 * d.rem(a, b); d.solve([\&](int i) { ans[i] =
 *   d.uf.find(qs[i].first) == d.uf.find(qs[i].second); });
 * Status: stress-tested against brute force
 */
#pragma once

#include "UnionFindRollback.h"

struct DynCon {
	RollbackUF uf;
	int q = 0;
	map<pii, vi> on; // start times of open copies
	vector<array<int, 4>> ed; // a, b, alive for queries [l, r)
	DynCon(int n) : uf(n) {}
	void add(int a, int b) { on[minmax(a, b)].push_back(q); }
	void rem(int a, int b) {
		vi& v = on[minmax(a, b)];
		ed.push_back({a, b, v.back(), q}); v.pop_back();
	}
	int query() { return q++; }
	template<class F> void rec(int l, int r, const vi& es, F& f) {
		int t = uf.time(), m = (l + r) / 2;
		vi sub;
		for (int i : es) {
			auto& e = ed[i];
			if (e[2] <= l && r <= e[3]) uf.join(e[0], e[1]);
			else if (e[2] < r && l < e[3]) sub.push_back(i);
		}
		if (r - l == 1) f(l);
		else rec(l, m, sub, f), rec(m, r, sub, f);
		uf.rollback(t);
	}
	template<class F> void solve(F f) {
		for (auto& [p, v] : on) for (int l : v)
			ed.push_back({p.first, p.second, l, q});
		vi es(sz(ed));
		iota(all(es), 0);
		if (q) rec(0, q, es, f);
	}
};
