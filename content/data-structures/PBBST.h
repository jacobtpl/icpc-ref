/**
 * Author: Siyong
 * Date: 2023-05-25
 * License: CC0
 * Source: me
 * Description: Persistent AVL tree, a sequential container
 *  with split/merge by 0-indexed position (range add, range
 *  sum). Versions work as in PersistentSegmentTree.h: nodes
 *  made since the last PAVL::step() are modified in place,
 *  older ones are copied first, so every tree that existed
 *  at a step() stays valid forever; split and + consume
 *  operands made after it. Nodes are never freed: call
 *  PAVL::reset() between testcases.
 *  All problem-specific code is in Node: apply puts a lazy
 *  tag on a subtree, up combines children, down pushes the
 *  tag (children are null if absent, and already copied),
 *  lazy says whether a tag is pending. Without lazies,
 *  delete add, apply and the body of down, and make lazy
 *  return 0.
 * Usage: PAVL a = PAVL(Node(1)) + PAVL(Node(2));
 *  PAVL::step(); // a is now a saved version
 *  auto [l, r] = a.split(1); // a intact
 *  r.get_root().apply(5); ll x = r.get_root().sum;
 *  PAVL b = l + r; // consumes l and r, not a
 * Time: O(log N) per operation
 * Status: stress-tested
 */
#pragma once

struct PAVL {
	struct Node {
		int t = T, s = 1, h = 1; array<int, 2> c{ -1, -1 };
		// Customize below
		ll val, sum, add = 0;
		Node(ll v = 0) : val(v), sum(v) {}
		void apply(ll x) { val += x, sum += x * s, add += x; }
		void up(const Node* a, const Node* b) {
			sum = val + (a ? a->sum : 0) + (b ? b->sum : 0);
		}
		bool lazy() const { return add != 0; }
		void down(Node* a, Node* b) {
			if (a) a->apply(add);
			if (b) b->apply(add);
			add = 0;
		}
	};
	static vector<Node> N;
	static int T;
	static int clone(int n) {
		if (n == -1 || N[n].t == T) return n;
		return N.push_back(N[n]), N.back().t = T, sz(N) - 1;
	}
	static Node* at(int n) { return n != -1 ? &N[n] : 0; }
	static int gh(int n) { return n != -1 ? N[n].h : 0; }
	static int gs(int n) { return n != -1 ? N[n].s : 0; }
	static void up(int n) {
		int a = N[n].c[0], b = N[n].c[1];
		N[n].h = max(gh(a), gh(b)) + 1;
		N[n].s = gs(a) + gs(b) + 1;
		N[n].up(at(a), at(b));
	}
	static int down(int n) {
		n = clone(n);
		if (N[n].lazy()) {
			int a = clone(N[n].c[0]), b = clone(N[n].c[1]);
			N[n].c = { a, b };
			N[n].down(at(a), at(b));
		}
		return n;
	}
	static int rotate(int n, int d) {
		n = clone(n); int o = down(N[n].c[d]);
		N[n].c[d] = N[o].c[!d], N[o].c[!d] = n;
		up(n), up(o);
		return o;
	}
	static int balance(int n) {
		assert(N[n].t == T); up(n);
		int diff = gh(N[n].c[0]) - gh(N[n].c[1]), d;
		if (diff >= 2) d = 0;
		else if (diff <= -2) d = 1;
		else return n;
		N[n].c[d] = down(N[n].c[d]);
		if (gh(N[N[n].c[d]].c[d]) + 1 < gh(N[n].c[d]))
			N[n].c[d] = rotate(N[n].c[d], !d);
		return rotate(n, d);
	}
	static int merge_root(int l, int n, int r) {
		if (gh(l) + 1 < gh(r))
			return r = down(r), N[r].c[0] = merge_root(l, n, N[r].c[0]), balance(r);
		else if (gh(r) + 1 < gh(l))
			return l = down(l), N[l].c[1] = merge_root(N[l].c[1], n, r), balance(l);
		else return N[n].c = { l, r }, balance(n);
	}
	static tuple<int, int> split(int n, int k) {
		if (n != -1) n = down(n);
		if (k == 0) return { -1, n };
		if (k == N[n].s) return { n, -1 };
		if (k <= gs(N[n].c[0])) {
			auto [l, r] = split(N[n].c[0], k);
			return { l, merge_root(r, n, N[n].c[1]) };
		} else {
			auto [l, r] = split(N[n].c[1], k - gs(N[n].c[0]) - 1);
			return { merge_root(N[n].c[0], n, l), r };
		}
	}
	static int merge(int l, int r) {
		if (r == -1) return l;
		auto [x, nr] = split(r, 1);
		return merge_root(l, clone(x), nr);
	}
	PAVL(int v) : root(v) {}
	int root;
	PAVL() : root(-1) {}
	PAVL(Node n) : root(sz(N)) { N.push_back(n); }
	friend PAVL operator+(PAVL a, PAVL b) { return merge(a.root, b.root); }
	tuple<PAVL, PAVL> split(int k) {
		auto [l, r] = split(root, k);
		return { PAVL(l), PAVL(r) };
	}
	static void step() { ++T; }
	static void reset() { N.clear(); T = 0; }
	Node& get_root() { return N[root = clone(root)]; } // nonempty only
};
typedef PAVL::Node Node;
vector<Node> PAVL::N;
int PAVL::T;
