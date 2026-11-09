#include "../utilities/template.h"
#include "../utilities/random.h"
#define pb push_back

#include "../../content/graph/MatroidIntersection.h"

// Independent brute-force oracles (not using the header's ok()).
bool isForest(int V, const vector<pii>& ed, const vi& S) {
	vi p(V);
	iota(all(p), 0);
	auto find = [&](int x) { while (p[x] != x) x = p[x]; return x; };
	for (int i : S) {
		int a = find(ed[i].first), b = find(ed[i].second);
		if (a == b) return false;
		p[a] = b;
	}
	return true;
}
bool colorOk(const vi& col, const vi& cap, const vi& S) {
	vi c = cap;
	for (int i : S) if (--c[col[i]] < 0) return false;
	return true;
}

// Uniform matroid of rank k: independent iff size <= k.
struct Uniform {
	int k, cur;
	void init(const vi& S) { cur = sz(S); }
	bool ok(int) { return cur < k; }
};
// Linear matroid over GF(2), vectors as bitmasks.
struct Linear {
	vi vec; vector<ll> basis;
	void init(const vi& S) {
		basis.assign(31, 0);
		for (int i : S) assert(ok(i)), add(vec[i]);
	}
	void add(ll x) {
		for (int b = 30; b >= 0; b--) if (x >> b & 1) {
			if (!basis[b]) { basis[b] = x; return; }
			x ^= basis[b];
		}
	}
	bool ok(int i) {
		ll x = vec[i];
		for (int b = 30; b >= 0; b--) if (x >> b & 1) {
			if (!basis[b]) return true;
			x ^= basis[b];
		}
		return false;
	}
};
bool linearOk(const vi& vec, const vi& S) {
	Linear l{vec, {}};
	l.basis.assign(31, 0);
	for (int i : S) {
		if (!l.ok(i)) return false;
		l.add(vec[i]);
	}
	return true;
}

template<class F1, class F2>
int bruteForce(int n, F1 ok1, F2 ok2) {
	int best = 0;
	rep(m,0,1 << n) {
		if (__builtin_popcount(m) <= best) continue;
		vi S;
		rep(i,0,n) if (m >> i & 1) S.pb(i);
		if (ok1(S) && ok2(S)) best = sz(S);
	}
	return best;
}

template<class A, class B, class F1, class F2>
void check(int n, A a, B b, F1 ok1, F2 ok2, int expected) {
	vi res = randBool() ? matroidIsect(n, a, b) : matroidIsect(n, b, a);
	assert(sz(res) == expected);
	rep(i,0,sz(res)) {
		assert(0 <= res[i] && res[i] < n);
		if (i) assert(res[i-1] < res[i]);
	}
	assert(ok1(res) && ok2(res));
}

struct Instance {
	int n, V;
	vector<pii> ed;
	vi col, cap, vec;
	int k;
};

Instance gen(int n) {
	Instance I;
	I.n = n;
	I.V = randRange(1, 8);
	rep(i,0,n) {
		int u = randRange(I.V), v = randRange(I.V);
		if (randRange(4)) v = randRange(I.V); else v = u; // some loops
		I.ed.pb({u, v});
	}
	int C = randRange(1, 6);
	I.cap.resize(C);
	rep(c,0,C) I.cap[c] = randRange(4);
	rep(i,0,n) I.col.pb(randRange(C));
	int dim = randRange(1, 6);
	rep(i,0,n) I.vec.pb(randRange(1 << dim)); // may be zero/duplicate
	I.k = randRange(n + 2);
	return I;
}

void testAllPairs(const Instance& I) {
	int n = I.n;
	auto gOk = [&](const vi& S) { return isForest(I.V, I.ed, S); };
	auto cOk = [&](const vi& S) { return colorOk(I.col, I.cap, S); };
	auto uOk = [&](const vi& S) { return sz(S) <= I.k; };
	auto lOk = [&](const vi& S) { return linearOk(I.vec, S); };
	Graphic g{I.V, I.ed};
	Colorful c{I.col, I.cap, {}};
	Uniform u{I.k, 0};
	Linear l{I.vec, {}};
	check(n, g, c, gOk, cOk, bruteForce(n, gOk, cOk));
	check(n, g, u, gOk, uOk, bruteForce(n, gOk, uOk));
	check(n, g, l, gOk, lOk, bruteForce(n, gOk, lOk));
	check(n, c, u, cOk, uOk, bruteForce(n, cOk, uOk));
	check(n, c, l, cOk, lOk, bruteForce(n, cOk, lOk));
	check(n, u, l, uOk, lOk, bruteForce(n, uOk, lOk));
	check(n, g, g, gOk, gOk, bruteForce(n, gOk, gOk));
	check(n, c, c, cOk, cOk, bruteForce(n, cOk, cOk));
	Instance J = gen(n);
	Graphic g2{J.V, J.ed};
	Colorful c2{J.col, J.cap, {}};
	auto gOk2 = [&](const vi& S) { return isForest(J.V, J.ed, S); };
	auto cOk2 = [&](const vi& S) { return colorOk(J.col, J.cap, S); };
	check(n, g, g2, gOk, gOk2, bruteForce(n, gOk, gOk2));
	check(n, c, c2, cOk, cOk2, bruteForce(n, cOk, cOk2));
}

void testEdgeCases() {
	// n = 0
	{
		Graphic g{3, {}};
		Colorful c{{}, {1, 2}, {}};
		Uniform u{5, 0};
		assert(matroidIsect(0, g, c).empty());
		assert(matroidIsect(0, c, g).empty());
		assert(matroidIsect(0, g, u).empty());
		Colorful c0{{}, {}, {}};
		assert(matroidIsect(0, c0, c0).empty());
	}
	// n = 1
	{
		Graphic g{2, {{0, 1}}};
		Graphic loop{1, {{0, 0}}};
		Colorful c1{{0}, {1}, {}};
		Colorful c0{{0}, {0}, {}};
		assert(matroidIsect(1, g, c1) == vi{0});
		assert(matroidIsect(1, c1, g) == vi{0});
		assert(matroidIsect(1, g, c0).empty());
		assert(matroidIsect(1, loop, c1).empty());
		assert(matroidIsect(1, c1, c1) == vi{0});
	}
	// all loops / all zero caps / zero-rank uniform
	{
		int n = 10;
		vector<pii> ed(n, {0, 0});
		Graphic g{1, ed};
		Colorful c{vi(n, 0), {5}, {}};
		assert(matroidIsect(n, g, c).empty());
		Graphic g2{n + 1, {}};
		rep(i,0,n) g2.ed.pb({i, i + 1});
		Colorful c0{vi(n, 0), {0}, {}};
		assert(matroidIsect(n, g2, c0).empty());
		Uniform u0{0, 0};
		assert(matroidIsect(n, g2, u0).empty());
		Colorful cbig{vi(n, 0), {100}, {}};
		vi all_(n); iota(all(all_), 0);
		assert(matroidIsect(n, g2, cbig) == all_);
	}
	// rainbow spanning tree: K_4 with one edge per colour -> 3
	{
		vector<pii> ed; vi col;
		rep(i,0,4) rep(j,i+1,4) ed.pb({i, j}), col.pb(sz(col));
		Graphic g{4, ed};
		Colorful c{col, vi(6, 1), {}};
		assert(sz(matroidIsect(6, g, c)) == 3);
		// all edges the same colour with cap 2 -> 2
		Colorful c2{vi(6, 0), {2}, {}};
		assert(sz(matroidIsect(6, g, c2)) == 2);
	}
	// large element indices / vertex ids (no overflow in int arithmetic)
	{
		int V = 1000000, n = 50;
		vector<pii> ed; vi col;
		rep(i,0,n) ed.pb({V - 1 - i, V - 2 - i}), col.pb(i % 7);
		Graphic g{V, ed};
		Colorful c{col, vi(7, INT_MAX), {}};
		assert(sz(matroidIsect(n, g, c)) == n);
	}
}

// Medium-size sanity: graphic x uniform = min(k, forest rank),
// colourful x uniform = min(k, sum min(cnt_c, cap_c)).
void testMedium() {
	rep(it,0,20) {
		int n = randRange(1, 400), V = randRange(1, 100);
		vector<pii> ed; vi col;
		rep(i,0,n) ed.pb({randRange(V), randRange(V)});
		int C = randRange(1, 20); vi cap(C);
		rep(c,0,C) cap[c] = randRange(5);
		rep(i,0,n) col.pb(randRange(C));
		int k = randRange(n + 1);
		Graphic g{V, ed};
		Colorful c{col, cap, {}};
		Uniform u{k, 0};
		vi full(n); iota(all(full), 0);
		int rankG = sz(matroidIsect(n, g, g));
		Graphic gg{V, ed}; gg.init({});
		int rk = 0; for (int i : full) if (gg.ok(i)) rk++, gg.uf.join(ed[i].first, ed[i].second);
		assert(rankG == rk);
		assert(sz(matroidIsect(n, g, u)) == min(k, rk));
		vi cnt(C); for (int x : col) cnt[x]++;
		int rc = 0; rep(i,0,C) rc += min(cnt[i], cap[i]);
		assert(sz(matroidIsect(n, c, u)) == min(k, rc));
		assert(sz(matroidIsect(n, c, c)) == rc);
		vi res = matroidIsect(n, g, c);
		assert(isForest(V, ed, res) && colorOk(col, cap, res));
		assert(sz(res) <= min(rk, rc));
	}
}

int main() {
	srand(2);
	testEdgeCases();
	rep(it,0,3000) testAllPairs(gen(randRange(0, 9)));
	rep(it,0,150) testAllPairs(gen(randRange(9, 13)));
	testMedium();
	cout << "Tests passed!" << endl;
}
