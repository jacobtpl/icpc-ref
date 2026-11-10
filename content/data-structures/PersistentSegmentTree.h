/**
 * Author: Jacob, Siyong
 * Date: 2026-11-08
 * License: CC0
 * Source: folklore
 * Description: Persistent lazy segment tree (range add, range sum) with
 * path copying. All problem-specific code lives in \texttt{Node}:
 * \texttt{Node()} must be the identity (empty segment), \texttt{new\_leaf}
 * makes a leaf from an input value, \texttt{empty\_range} is the untouched
 * segment $[lo, hi)$ (nodes are created lazily from it), \texttt{apply}
 * applies a lazy tag to the segment $[lo, hi)$, \texttt{up} combines
 * children, \texttt{down} pushes the tag to child copies.
 * Bounds are $[lo, hi)$ and passed as arguments; root 0 is the untouched
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
	static Node new_leaf(ll val) { Node n; n.v = val; return n; }
	// empty_range(lo, hi) = untouched segment; must equal
	// up(empty_range(lo, m), empty_range(m, hi), lo, hi)
	static Node empty_range(int, int) { return Node(); }
	void apply(ll x, int lo, int hi) {
		v += x * (hi - lo), lz += x;
	}
	void up(const Node& a, const Node& b, int lo, int hi) {
		v = a.v + b.v;
	}
	void down(Node& a, Node& b, int lo, int hi) {
		int m = lo + (hi - lo) / 2;
		a.apply(lz, lo, m), b.apply(lz, m, hi), lz = 0;
	}
};
vector<Node> N(1); int T = 0;
void step() { ++T; }
void reset() { N.assign(1, Node()); T = 0; }
int make(Node n) { n.t = T; N.push_back(n); return sz(N) - 1; }
Node get(int x, int lo, int hi) {
	return x ? N[x] : Node::empty_range(lo, hi);
}
int clone(int x, int lo, int hi) {
	if (!x) return make(Node::empty_range(lo, hi));
	return N[x].t == T ? x : make(N[x]);
}
void up(int x, int lo, int hi) {
	int m = lo + (hi - lo) / 2;
	N[x].up(get(N[x].c[0], lo, m), get(N[x].c[1], m, hi), lo, hi);
}
void down(int x, int lo, int hi) { // x must be cloned
	int m = lo + (hi - lo) / 2;
	int a = clone(N[x].c[0], lo, m), b = clone(N[x].c[1], m, hi);
	N[x].c = {a, b};
	N[x].down(N[a], N[b], lo, hi);
}
template<class T> int build(vector<T>& v, int lo, int hi) {
	if (lo + 1 == hi) return make(Node::new_leaf(v[lo]));
	int x = make(Node()), m = lo + (hi - lo) / 2;
	N[x].c = {build(v, lo, m), build(v, m, hi)};
	return up(x, lo, hi), x;
}
int add(int x, int lo, int hi, int l, int r, ll val) {
	if (r <= lo || hi <= l) return x;
	x = clone(x, lo, hi);
	if (l <= lo && hi <= r) return N[x].apply(val, lo, hi), x;
	down(x, lo, hi);
	int m = lo + (hi - lo) / 2;
	N[x].c[0] = add(N[x].c[0], lo, m, l, r, val);
	N[x].c[1] = add(N[x].c[1], m, hi, l, r, val);
	return up(x, lo, hi), x;
}
Node query(Node x, int lo, int hi, int l, int r) {
	if (r <= lo || hi <= l) return Node();
	if (l <= lo && hi <= r) return x;
	int m = lo + (hi - lo) / 2;
	Node a = get(x.c[0], lo, m), b = get(x.c[1], m, hi), res;
	x.down(a, b, lo, hi);
	res.up(query(a, lo, m, l, r), query(b, m, hi, l, r), lo, hi);
	return res;
}
Node query(int x, int lo, int hi, int l, int r) {
	return query(get(x, lo, hi), lo, hi, l, r);
}
}
