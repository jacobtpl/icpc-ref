#include "../utilities/template.h"

#include "../../content/graph/GeneralMatching.h"

mt19937 rng(777);
int ri(int a, int b) { return (int)(rng() % (unsigned)(b - a + 1)) + a; }

// Maximum matching size by bitmask DP.
int brute(int n, const vector<vi>& g) {
	vi dp(1 << n);
	rep(m,1,1<<n) {
		int i = __builtin_ctz(m);
		dp[m] = dp[m ^ (1 << i)];
		rep(j,i+1,n) if ((m >> j & 1) && g[i][j])
			dp[m] = max(dp[m], dp[m ^ (1 << i) ^ (1 << j)] + 1);
	}
	return dp[(1 << n) - 1];
}

void check(int n, vector<pii> ed) {
	vector<vi> g(n, vi(n));
	for (pii e : ed) g[e.first][e.second] = g[e.second][e.first] = 1;
	auto r = generalMatching(n, ed);
	vi used(n);
	for (pii e : r) {
		assert(0 <= e.first && e.first < n && 0 <= e.second && e.second < n);
		assert(g[e.first][e.second]);
		assert(!used[e.first] && !used[e.second]);
		used[e.first] = used[e.second] = 1;
	}
	assert(sz(r) == brute(n, g));
}

void test(int maxn, int iters) {
	rep(it,0,iters) {
		int n = ri(0, maxn), m = n ? ri(0, n * n) : 0, dens = ri(0, 100);
		vector<pii> ed;
		rep(e,0,m) { // duplicates and both orientations allowed
			int a = ri(0, n-1), b = ri(0, n-1);
			if (a != b && ri(0, 99) < dens) ed.push_back({a, b});
		}
		check(n, ed);
	}
}

int main() {
	srand(3);
	check(0, {}); check(1, {}); check(2, {}); check(2, {{1, 0}});
	check(3, {{0, 1}, {1, 2}, {2, 0}});
	rep(n,1,13) { // cliques, cycles, paths, stars
		vector<pii> cl, cy, pa, st;
		rep(i,0,n) rep(j,0,i) cl.push_back({i, j});
		rep(i,0,n) if (n > 2) cy.push_back({i, (i + 1) % n});
		rep(i,1,n) pa.push_back({i - 1, i}), st.push_back({0, i});
		check(n, cl); check(n, cy); check(n, pa); check(n, st);
	}
	test(4, 30000);
	test(8, 20000);
	test(12, 10000);
	test(16, 300);
#ifdef SELF_LOOPS // undocumented precondition: a self-loop trips assert(r % 2 == 0)
	check(1, {{0, 0}});
#endif
	cout<<"Tests passed!"<<endl;
}
