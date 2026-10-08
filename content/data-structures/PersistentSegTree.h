/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: folklore (path copying); written from scratch
 * Description: Persistent sum segment tree over positions $[0, n)$,
 *  $1 \le n \le 2^{30}$, all initially 0. A version is a root index;
 *  root 0 is the all-zero tree. \texttt{upd} adds $x$ at position $p$
 *  and returns the new root, old roots stay valid. \texttt{query} sums
 *  $[b, e)$ in any version. \texttt{kth(a, b, k)} is the least $p$ with
 *  more than $k$ in $[0, p]$ of version $b$ minus version $a$; needs
 *  all point differences $\ge 0$ and $0 \le k <$ total difference.
 *  The tree is implicit, so $n$ can be the value range. Each update
 *  adds $\le \lceil \log_2 n \rceil + 1$ nodes (16 bytes each): reserve.
 * Usage:
 *  PST t(m); vi r{0}; // values a[i] in [0, m)
 *  for (int x : a) r.push_back(t.upd(r.back(), x, 1));
 *  t.kth(r[lo], r[hi], k); // k-th (0-indexed) smallest in a[lo, hi)
 *  t.query(r[hi], x, y) - t.query(r[lo], x, y); // count in [x, y)
 * Time: $O(\log n)$ per operation.
 * Memory: $O(\log n)$ per update.
 * Status: stress-tested
 */
#pragma once

struct PST {
	typedef ll T;
	struct Node { int l, r; T v; };
	vector<Node> t; int n;
	PST(int n) : t(1), n(n) {}
	int upd(int o, int p, T x) { // a[p] += x in version o
		int res = sz(t), lo = 0, hi = n;
		for (;;) {
			Node c = t[o]; c.v += x;
			if (hi - lo == 1) { t.push_back(c); break; }
			int m = (lo + hi) / 2;
			if (p < m) o = c.l, c.l = sz(t) + 1, hi = m;
			else o = c.r, c.r = sz(t) + 1, lo = m;
			t.push_back(c);
		}
		return res;
	}
	T q(int o, int b, int e, int lo, int hi) {
		if (!o || b >= hi || e <= lo) return 0;
		if (b <= lo && hi <= e) return t[o].v;
		int m = (lo + hi) / 2;
		return q(t[o].l, b, e, lo, m) + q(t[o].r, b, e, m, hi);
	}
	T query(int o, int b, int e) { return q(o, b, e, 0, n); }
	int kth(int a, int b, T k) {
		int lo = 0, hi = n;
		while (hi - lo > 1) {
			int m = (lo + hi) / 2;
			T c = t[t[b].l].v - t[t[a].l].v;
			if (k < c) a = t[a].l, b = t[b].l, hi = m;
			else k -= c, a = t[a].r, b = t[b].r, lo = m;
		}
		return lo;
	}
};
