#include "../utilities/template.h"

#include "../../content/graph/DFSMatching.h"

mt19937 rng(12345);
int rnd(int lo, int hi) { return (int)(rng() % (unsigned)(hi - lo + 1)) + lo; }

// O(n * 2^m) oracle: maximum matching by DP over subsets of the right side.
int brute(vector<vi>& g, int m) {
	vi dp(1 << m, -1), nd;
	dp[0] = 0;
	for (auto& adj : g) {
		nd = dp;
		rep(mask,0,1 << m) if (dp[mask] >= 0)
			for (int j : adj) if (!(mask >> j & 1))
				nd[mask | 1 << j] = max(nd[mask | 1 << j], dp[mask] + 1);
		dp = nd;
	}
	return *max_element(all(dp));
}

void check(vector<vi>& g, int m, int expected) {
	vi btoa(m, -1);
	int res = dfsMatching(g, btoa);
	assert(res == expected);
	int cnt = 0;
	vi usedL(sz(g));
	rep(j,0,m) if (btoa[j] != -1) {
		int a = btoa[j];
		assert(0 <= a && a < sz(g));
		assert(!usedL[a]++);
		assert(count(all(g[a]), j));
		cnt++;
	}
	assert(cnt == res);
}

int main() {
	// edge cases: empty sides, no edges
	rep(n,0,4) rep(m,0,4) {
		vector<vi> g(n);
		check(g, m, 0);
	}
	// small random, with duplicate edges and isolated vertices
	rep(it,0,60000) {
		int n = rnd(1, 8), m = rnd(1, 8);
		int e = rnd(0, n * m + 3);
		if (it % 3 == 0) e = rnd(0, max(n, m));
		vector<vi> g(n);
		rep(i,0,e) g[rnd(0, n-1)].push_back(rnd(0, m-1));
		check(g, m, brute(g, m));
	}
	// complete bipartite
	rep(n,1,9) rep(m,1,9) {
		vector<vi> g(n);
		rep(i,0,n) rep(j,0,m) g[i].push_back(j);
		check(g, m, min(n, m));
	}
	// medium random with a planted perfect matching
	rep(it,0,200) {
		int n = rnd(50, 300);
		vi perm(n);
		iota(all(perm), 0);
		shuffle(all(perm), rng);
		vector<vi> g(n);
		rep(i,0,n) {
			g[i].push_back(perm[i]);
			rep(k,0,rnd(0, 3)) g[i].push_back(rnd(0, n-1));
			shuffle(all(g[i]), rng);
		}
		check(g, n, n);
	}
	// one augmenting path through every vertex (recursion depth n)
	{
		int n = 20000;
		vector<vi> g(n);
		rep(i,0,n-1) g[i] = {i + 1, i};
		g[n-1] = {n - 1};
		check(g, n, n);
	}
#ifdef DEEP_RECURSION
	// Opt-in: find() recurses once per vertex of the augmenting path, which
	// needs more than the default 8 MB stack here (run with ulimit -s unlimited).
	{
		int n = 150000;
		vector<vi> g(n);
		rep(i,0,n-1) g[i] = {i + 1, i};
		g[n-1] = {n - 1};
		check(g, n, n);
	}
#endif
	cout << "Tests passed!" << endl;
}
