#include "../utilities/template.h"

#include "../../content/graph/hopcroftKarp.h"

// Oracle: plain Kuhn augmenting paths.
bool kuhn(int a, vector<vi>& g, vi& mt, vi& vis) {
	for (int b : g[a]) if (!vis[b]) {
		vis[b] = 1;
		if (mt[b] == -1 || kuhn(mt[b], g, mt, vis)) return mt[b] = a, 1;
	}
	return 0;
}
int slow(vector<vi>& g, int m) {
	vi mt(m, -1), vis; int r = 0;
	rep(a,0,sz(g)) vis.assign(m, 0), r += kuhn(a, g, mt, vis);
	return r;
}
// Oracle 2: bitmask DP over the right side (m <= 12).
int dp(vector<vi>& g, int m) {
	vi best(1 << m, -1); best[0] = 0;
	for (auto& li : g) {
		vi nb = best;
		rep(s,0,1 << m) if (best[s] >= 0) for (int b : li) if (!(s >> b & 1))
			nb[s | 1 << b] = max(nb[s | 1 << b], best[s] + 1);
		best = nb;
	}
	return *max_element(all(best));
}

mt19937 rng(2024);
int rnd(int lo, int hi) { return lo + (int)(rng() % (unsigned)(hi - lo + 1)); }

void check(vector<vi>& g, int m, int expect) {
	vi btoa(m, -1);
	int res = hopcroftKarp(g, btoa);
	assert(res == expect);
	int cnt = 0; vi usedA(sz(g));
	rep(b,0,m) if (btoa[b] != -1) {
		int a = btoa[b]; cnt++;
		assert(0 <= a && a < sz(g) && !usedA[a]++);
		assert(find(all(g[a]), b) != g[a].end());
	}
	assert(cnt == res);
}

int main() {
	{ vector<vi> g; check(g, 0, 0); check(g, 3, 0); }
	{ vector<vi> g(3); check(g, 0, 0); check(g, 2, 0); }
	{ vector<vi> g{{0, 0, 0}}; check(g, 1, 1); }
	rep(it,0,200000) {
		int n = rnd(0, 7), m = rnd(0, 7), p = rnd(0, 100);
		vector<vi> g(n);
		rep(i,0,n) rep(j,0,m) rep(k,0,2) if (rnd(0, 99) < (k ? p / 4 : p)) g[i].push_back(j); // multi-edges
		rep(i,0,n) shuffle(all(g[i]), rng);
		int a = slow(g, m);
		assert(a == dp(g, m));
		check(g, m, a);
	}
	rep(it,0,20000) {
		int n = rnd(1, 60), m = rnd(1, 60), e = rnd(0, 3 * (n + m));
		vector<vi> g(n);
		rep(i,0,e) g[rnd(0, n - 1)].push_back(rnd(0, m - 1));
		check(g, m, slow(g, m));
	}
	// long augmenting paths: a_i - {b_i, b_{i+1}}, listed so that greedy choices are bad
	rep(it,0,200) {
		int n = rnd(1, 300);
		vector<vi> g(n);
		rep(i,0,n) { g[i] = {i + 1, i}; if (it % 2) swap(g[i][0], g[i][1]); }
		g.push_back({0});
		if (it % 4 < 2) rotate(g.begin(), g.end() - 1, g.end());
		check(g, n + 1, n + 1);
	}
	// btoa may be pre-seeded with a valid partial matching
	rep(it,0,20000) {
		int n = rnd(1, 8), m = rnd(1, 8);
		vector<vi> g(n);
		rep(i,0,n) rep(j,0,m) if (rnd(0, 2) == 0) g[i].push_back(j);
		vi btoa(m, -1), used(n);
		rep(j,0,m) { int i = rnd(0, n - 1); if (!used[i] && count(all(g[i]), j) && rnd(0, 1)) used[i] = 1, btoa[j] = i; }
		int pre = m - (int)count(all(btoa), -1);
		int res = hopcroftKarp(g, btoa);
		assert(pre + res == slow(g, m));
		assert(m - (int)count(all(btoa), -1) == pre + res);
	}
	cout << "Tests passed!" << endl;
}
