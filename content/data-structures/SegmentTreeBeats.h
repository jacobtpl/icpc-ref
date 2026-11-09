/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: Ji Yuhao's 2016 "Segment Tree Beats" (via the Codeforces blog by jiry\_2); own implementation
 * Description: Segment tree beats on a 0-indexed array, half-open ranges $[L, R)$.
 * \texttt{update(L,R,v,k)}: $k=0$ sets $a_i = \min(a_i, v)$, $k=1$ sets $a_i = \max(a_i, v)$, $k=2$ adds $v$.
 * \texttt{query(L,R,k)}: $k=0$ max, $k=1$ min, $k=2$ sum; an empty range gives $-\infty$, $\infty$, $0$.
 * Each node keeps the max, strict second max and count of max in \texttt{m[0]}, \texttt{m2[0]}, \texttt{c[0]},
 * and the same for the negated min in index 1.
 * Range sums must fit in ll, and $\max|a_i| + \sum|v_{add}| < 2^{62}$ (also for $v$ in chmin/chmax):
 * lazy adds pile up below a clamped node, so node sums are unsigned and may wrap until pushed.
 * Time: $O(\log N)$ per query, amortised $O(\log^2 N)$ per update ($O(\log N)$ if there are no adds).
 * Usage: Beats t(a); t.update(l, r, 5, 0); t.query(l, r, 2);
 * Status: stress-tested
 */
#pragma once

struct Beats {
	static constexpr ll inf = LLONG_MAX;
	struct N { uint64_t s; ll ad, m[2], m2[2]; int c[2]; };
	int n; vector<N> t;
	Beats(vector<ll>& a) : n(sz(a)), t(4 * n) {
		if (n) build(a, 1, 0, n);
	}
	void build(vector<ll>& a, int x, int lo, int hi) {
		if (lo + 1 == hi) {
			ll v = a[lo];
			t[x] = {(uint64_t)v, 0, {v, -v}, {-inf, -inf}, {1, 1}};
			return;
		}
		int mid = (lo + hi) / 2;
		build(a, 2*x, lo, mid), build(a, 2*x+1, mid, hi);
		pull(x);
	}
	void pull(int x) {
		N &a = t[2*x], &b = t[2*x+1], &c = t[x];
		c.s = a.s + b.s;
		rep(k,0,2) {
			c.m[k] = max(a.m[k], b.m[k]);
			c.m2[k] = max(a.m2[k], b.m2[k]), c.c[k] = 0;
			for (N* p : {&a, &b})
				if (p->m[k] == c.m[k]) c.c[k] += p->c[k];
				else c.m2[k] = max(c.m2[k], p->m[k]);
		}
	}
	void add(int x, ll v, int len) {
		N& a = t[x]; a.s += (uint64_t)v * len, a.ad += v;
		rep(k,0,2) {
			a.m[k] += k ? -v : v;
			if (a.m2[k] != -inf) a.m2[k] += k ? -v : v;
		}
	}
	void cut(int x, int k, ll v) { // m2[k] < v < m[k]
		N& a = t[x];
		a.s += uint64_t(v - a.m[k]) * a.c[k] * (k ? -1 : 1);
		if (a.m[!k] == -a.m[k]) a.m[!k] = -v;
		else if (a.m2[!k] == -a.m[k]) a.m2[!k] = -v;
		a.m[k] = v;
	}
	void push(int x, int lo, int hi) {
		int mid = (lo + hi) / 2;
		rep(i,0,2) {
			int y = 2*x + i;
			add(y, t[x].ad, i ? hi - mid : mid - lo);
			rep(k,0,2) if (t[y].m[k] > t[x].m[k])
				cut(y, k, t[x].m[k]);
		}
		t[x].ad = 0;
	}
	void upd(int L, int R, ll v, int k, int x, int lo, int hi) {
		if (R <= lo || hi <= L || (k < 2 && t[x].m[k] <= v))
			return;
		if (L <= lo && hi <= R && (k == 2 || t[x].m2[k] < v))
			return k == 2 ? add(x, v, hi - lo) : cut(x, k, v);
		int mid = (lo + hi) / 2; push(x, lo, hi);
		upd(L, R, v, k, 2*x, lo, mid);
		upd(L, R, v, k, 2*x+1, mid, hi), pull(x);
	}
	ll qry(int L, int R, int k, int x, int lo, int hi) {
		if (R <= lo || hi <= L) return k == 2 ? 0 : -inf;
		if (L <= lo && hi <= R)
			return k == 2 ? (ll)t[x].s : t[x].m[k];
		int mid = (lo + hi) / 2; push(x, lo, hi);
		ll a = qry(L, R, k, 2*x, lo, mid);
		ll b = qry(L, R, k, 2*x+1, mid, hi);
		return k == 2 ? a + b : max(a, b);
	}
	void update(int L, int R, ll v, int k) {
		upd(L, R, k == 1 ? -v : v, k, 1, 0, n);
	}
	ll query(int L, int R, int k) {
		ll r = qry(L, R, k, 1, 0, n);
		return k == 1 ? -r : r;
	}
};
