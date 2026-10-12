/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: own work; idea from
 *  https://cp-algorithms.com/geometry/convex_hull_trick.html
 * Description: Li Chao tree over integer $x \in [lo, hi)$ with
 *  nodes created on demand. Add lines $kx+m$, optionally only
 *  on $x \in [a, b)$ (clipped to $[lo, hi)$), and query the
 *  minimum at a point; \texttt{inf} if no line covers it.
 *  For maximum, negate $k$, $m$ and the answer.
 *  $hi - lo$ and $kx+m$ for every $x$ the line is added on
 *  must fit in ll. Lines may be added in any order.
 *  Memory: one node per line, $O(\log C)$ nodes per segment
 *  ($\approx 40$ at $C = 2^{31}$), so keep $[lo, hi)$ tight.
 * Usage: LiChao t(-1e9, 1e9 + 1); t.add(k, m);
 *  t.add(k, m, a, b); t.query(x);
 * Time: O(\log C) per line and query, O(\log^2 C) per segment,
 *  $C = hi - lo$
 * Status: stress-tested
 */
#pragma once

struct LiChao {
	static const ll inf = LLONG_MAX;
	struct Node {
		ll k, m; Node *l = 0, *r = 0;
		ll f(ll x) { return k * x + m; }
	} *root = 0;
	ll lo, hi;
	LiChao(ll lo, ll hi) : lo(lo), hi(hi) {}
	void add(Node*& n, ll l, ll r, ll a, ll b, ll k, ll m) {
		if (b <= l || r <= a) return;
		ll mid = l + (r - l) / 2;
		if (a > l || r > b) {
			if (!n) n = new Node{0, inf};
			add(n->l, l, mid, a, b, k, m);
			add(n->r, mid, r, a, b, k, m);
		} else if (!n) n = new Node{k, m};
		else {
			bool L = k * l + m < n->f(l);
			bool M = k * mid + m < n->f(mid);
			if (M) swap(k, n->k), swap(m, n->m);
			if (r - l == 1) return;
			if (L != M) add(n->l, l, mid, a, b, k, m);
			else add(n->r, mid, r, a, b, k, m);
		}
	}
	void add(ll k, ll m) { add(root, lo, hi, lo, hi, k, m); }
	void add(ll k, ll m, ll a, ll b) { // x in [a, b)
		add(root, lo, hi, a, b, k, m);
	}
	ll query(ll x) { // lo <= x < hi
		ll l = lo, r = hi, res = inf;
		for (Node* n = root; n;) {
			res = min(res, n->f(x));
			ll mid = l + (r - l) / 2;
			if (x < mid) n = n->l, r = mid;
			else n = n->r, l = mid;
		}
		return res;
	}
};
