/**
 * Author: Devin
 * Date: 2026-11-17
 * License: CC0
 * Source: folklore (weighted / potential DSU); UCup 5 Stage 6 F
 * Description: DSU where every node has an unknown value p(x) and each join
 *  states a difference: join(a, b, v) asserts p(a) - p(b) = v. Returns 1 if
 *  two sets were merged, 0 if the statement was already implied, -1 if it
 *  contradicts earlier ones (nothing changes). diff(a, b) = p(a) - p(b),
 *  valid only if same(a, b).
 *  Group: change T, add and neg. Xor (parity / path xor): add = a ^ b,
 *  neg = a. Bipartiteness: xor with v = 1.
 *  No path compression, so rollback works; drop st/time/rollback if unused.
 * Time: O(log N) per operation
 * Usage: int t = uf.time(); uf.join(a, b, 5); uf.rollback(t);
 * Status: stress-tested against brute force (sum and xor, with rollback)
 */
#pragma once

struct PotentialUF {
	typedef ll T;
	static T add(T a, T b) { return a + b; }
	static T neg(T a) { return -a; }

	vi e; vector<T> w; vector<pii> st; // w[x] = p(x) - p(parent)
	PotentialUF(int n) : e(n, -1), w(n) {}
	pair<int, T> find(int x) { // {root, p(x) - p(root)}
		T d{};
		for (; e[x] >= 0; x = e[x]) d = add(d, w[x]);
		return {x, d};
	}
	bool same(int a, int b) { return find(a).first == find(b).first; }
	int size(int x) { return -e[find(x).first]; }
	T diff(int a, int b) {
		return add(find(a).second, neg(find(b).second));
	}
	int join(int a, int b, T v) {
		auto [ra, da] = find(a); auto [rb, db] = find(b);
		if (ra == rb) return add(da, neg(db)) == v ? 0 : -1;
		T x = add(add(v, db), neg(da)); // p(ra) - p(rb)
		if (e[ra] < e[rb]) swap(ra, rb), x = neg(x);
		st.push_back({ra, e[ra]});
		e[rb] += e[ra]; e[ra] = rb; w[ra] = x;
		return 1;
	}
	int time() { return sz(st); }
	void rollback(int t) {
		for (; sz(st) > t; st.pop_back()) {
			auto [c, old] = st.back();
			e[e[c]] -= old; e[c] = old;
		}
	}
};
