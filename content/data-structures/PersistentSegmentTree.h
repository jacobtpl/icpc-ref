/**
 * Author: Jacob, Siyong
 * Date: 2026-11-08
 * License: CC0
 * Source: folklore
 * Description: Persistent lazy segment tree (range add, range sum) with
 * path copying. All problem-specific code lives in \texttt{Node}:
 * \texttt{Node()} must be the identity (empty segment), \texttt{apply}
 * applies a lazy tag to a segment of length \texttt{len}, \texttt{up}
 * combines children, \texttt{down} pushes the tag to child copies.
 * Bounds are $[lo, hi)$ and passed as arguments; root 0 is the all-identity
 * tree, so \texttt{build} is optional and coordinates up to $10^9$ work
 * without compression. Mutating ops return the new root. Nodes created
 * since the last \texttt{step()} are modified in place (the returned root
 * may equal the input); older versions are never changed. Nodes are never
 * freed: \texttt{reset()} between testcases. For point updates only,
 * delete \texttt{lz}, \texttt{down} and the calls to \texttt{down}.
 * Usage: int r0 = pst::build(v, 0, n); // or r0 = 0
 *  pst::step(); int r1 = pst::add(r0, 0, n, l, r, 5); // r0 intact
 *  ll s = pst::query(r1, 0, n, l, r).v;
 * Time: O(\log N) per op; $\le 4\log N$ new nodes per range op,
 * $2 \log N$ per point op (half that without lazy).
 * Status: stress-tested
 */
#pragma once

namespace pst {
struct Node {
	array<int, 2> c = {0, 0}; int t = -1;
	ll v = 0, lz = 0; /* Customize below */
	void apply(ll x, int len) { v += x * len, lz += x; }
	void up(const Node& a, const Node& b) { v = a.v + b.v; }
	void down(Node& a, Node& b, int la, int lb) {
		a.apply(lz, la), b.apply(lz, lb), lz = 0;
	}
};
vector<Node> N(1); int T = 0;
void step() { ++T; }
void reset() { N.assign(1, Node()); T = 0; }
int clone(int x) {
	if (N[x].t == T) return x;
	N.push_back(N[x]); N.back().t = T;
	return sz(N) - 1;
}
void up(int x) { N[x].up(N[N[x].c[0]], N[N[x].c[1]]); }
void down(int x, int lo, int hi) { // x must be cloned
	int m = (lo + hi) / 2;
	int a = clone(N[x].c[0]), b = clone(N[x].c[1]);
	N[x].c = {a, b};
	N[x].down(N[a], N[b], m - lo, hi - m);
}
int build(vector<ll>& v, int lo, int hi) {
	int x = clone(0), m = (lo + hi) / 2;
	if (lo + 1 == hi) return N[x].v = v[lo], x;
	N[x].c = {build(v, lo, m), build(v, m, hi)};
	return up(x), x;
}
int add(int x, int lo, int hi, int l, int r, ll val) {
	if (r <= lo || hi <= l) return x;
	x = clone(x);
	if (l <= lo && hi <= r) return N[x].apply(val, hi - lo), x;
	down(x, lo, hi);
	int m = (lo + hi) / 2;
	N[x].c[0] = add(N[x].c[0], lo, m, l, r, val);
	N[x].c[1] = add(N[x].c[1], m, hi, l, r, val);
	return up(x), x;
}
Node query(Node x, int lo, int hi, int l, int r) {
	if (r <= lo || hi <= l) return Node();
	if (l <= lo && hi <= r) return x;
	int m = (lo + hi) / 2;
	Node a = N[x.c[0]], b = N[x.c[1]], res;
	x.down(a, b, m - lo, hi - m);
	res.up(query(a, lo, m, l, r), query(b, m, hi, l, r));
	return res;
}
Node query(int x, int lo, int hi, int l, int r) {
	return query(N[x], lo, hi, l, r);
}
}
