#include "../utilities/template.h"

#include "../../content/data-structures/PersistentSegmentTree.h"

// Copy of the header (#pragma once prevents re-including it) with
// the untouched array being a[i] = i, to exercise empty_range.
namespace pos { namespace pst {
struct Node {
	array<int, 2> c = {0, 0}; int t = -1;
	ll v = 0, lz = 0;
	static Node empty_range(int lo, int hi) {
		Node n; n.v = (ll)(lo + hi - 1) * (hi - lo) / 2; return n;
	}
	void apply(ll x, int lo, int hi) {
		v += x * (hi - lo), lz += x;
	}
	void up(const Node& a, const Node& b) { v = a.v + b.v; }
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
	N[x].up(get(N[x].c[0], lo, m), get(N[x].c[1], m, hi));
}
void down(int x, int lo, int hi) { // x must be cloned
	int m = lo + (hi - lo) / 2;
	int a = clone(N[x].c[0], lo, m), b = clone(N[x].c[1], m, hi);
	N[x].c = {a, b};
	N[x].down(N[a], N[b], lo, hi);
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
	res.up(query(a, lo, m, l, r), query(b, m, hi, l, r));
	return res;
}
Node query(int x, int lo, int hi, int l, int r) {
	return query(get(x, lo, hi), lo, hi, l, r);
}
}}

static unsigned R = 1231231;
int ra() {
	R *= 574841;
	R += 14;
	return (int)(R >> 1);
}

ll brute(vector<ll>& v, int l, int r) {
	ll s = 0;
	rep(k,l,r) s += v[k];
	return s;
}

void check(int root, vector<ll>& v, int n, int l, int r) {
	assert(pst::query(root, 0, n, l, r).v == brute(v, l, r));
}

void run(int n, bool useBuild) {
	pst::reset();
	vector<ll> init(n);
	if (useBuild) for (ll& x : init) x = ra() % 21 - 10;
	vector<int> init32(all(init));
	int r0 = useBuild ? pst::build(init32, 0, n) : 0;
	vector<pair<int, vector<ll>>> vers = {{r0, init}};
	rep(it,0,3000) {
		pst::step();
		auto [root, v] = vers[ra() % sz(vers)];
		int ops = ra() % 4;
		rep(o,0,ops) {
			int l = ra() % (n + 1), r = ra() % (n + 1);
			if (l > r) swap(l, r);
			ll x = ra() % 21 - 10;
			root = pst::add(root, 0, n, l, r, x);
			rep(k,l,r) v[k] += x;
			int ql = ra() % (n + 1), qr = ra() % (n + 1);
			if (ql > qr) swap(ql, qr);
			check(root, v, n, ql, qr);
		}
		vers.emplace_back(root, v);
		// older versions must be intact
		rep(k,0,3) {
			auto& [oroot, ov] = vers[ra() % sz(vers)];
			int ql = ra() % (n + 1), qr = ra() % (n + 1);
			if (ql > qr) swap(ql, qr);
			check(oroot, ov, n, ql, qr);
		}
	}
	for (auto& [oroot, ov] : vers) check(oroot, ov, n, 0, n);
}

void sparse() {
	// coordinates up to 1e9 without build
	pst::reset();
	const int C = 1000000000;
	map<int, ll> v;
	vector<int> roots = {0};
	vector<map<int, ll>> vs = {v};
	rep(it,0,2000) {
		pst::step();
		int i = ra() % sz(roots);
		int p = ra() % 2 ? ra() % C : ra() % 50;
		ll x = ra() % 21 - 10;
		int root = pst::add(roots[i], 0, C, p, p + 1, x);
		auto nv = vs[i];
		nv[p] += x;
		roots.push_back(root), vs.push_back(nv);
		int j = ra() % sz(roots);
		int l = ra() % C, r = ra() % C;
		if (l > r) swap(l, r);
		ll s = 0;
		for (auto& [q, y] : vs[j]) if (l <= q && q < r) s += y;
		assert(pst::query(roots[j], 0, C, l, r).v == s);
		assert(pst::query(roots[j], 0, C, 0, C).v
			== pst::query(roots[j], 0, C, 0, 50).v
			+ pst::query(roots[j], 0, C, 50, C).v);
	}
}

void lazyInit(int n) {
	// untouched a[i] = i, created lazily via empty_range
	namespace pst = pos::pst;
	pst::reset();
	vector<ll> init(n);
	iota(all(init), 0);
	vector<pair<int, vector<ll>>> vers = {{0, init}};
	rep(it,0,2000) {
		pst::step();
		auto [root, v] = vers[ra() % sz(vers)];
		rep(o,0,ra() % 3) {
			int l = ra() % (n + 1), r = ra() % (n + 1);
			if (l > r) swap(l, r);
			ll x = ra() % 21 - 10;
			root = pst::add(root, 0, n, l, r, x);
			rep(k,l,r) v[k] += x;
		}
		vers.emplace_back(root, v);
		auto& [oroot, ov] = vers[ra() % sz(vers)];
		int ql = ra() % (n + 1), qr = ra() % (n + 1);
		if (ql > qr) swap(ql, qr);
		assert(pst::query(oroot, 0, n, ql, qr).v == brute(ov, ql, qr));
	}
	// huge coordinates, nothing materialized
	const int C = 1000000000;
	int l = ra() % C, r = ra() % C;
	if (l > r) swap(l, r);
	assert(pst::query(0, 0, C, l, r).v == (ll)(l + r - 1) * (r - l) / 2);
	int root = pst::add(0, 0, C, l, r, 3);
	assert(pst::query(root, 0, C, 0, C).v
		== (ll)(C - 1) * C / 2 + 3LL * (r - l));
}

int main() {
	rep(n,1,12) run(n, true), run(n, false);
	run(50, true);
	run(200, false);
	sparse();
	rep(n,1,12) lazyInit(n);
	lazyInit(100);
	cout << "Tests passed!" << endl;
}
