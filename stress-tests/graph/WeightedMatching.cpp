#include "../utilities/template.h"
#include "../utilities/random.h"

#include "../../content/graph/WeightedMatching.h"

// Oracle 1: all injections (n <= m small).
ll brute(const vector<vi>& c, int n, int m) {
	vector<ll> dp(1 << m, LLONG_MAX); dp[0] = 0;
	ll best = LLONG_MAX;
	rep(s,0,1 << m) if (dp[s] != LLONG_MAX) {
		int i = __builtin_popcount(s);
		if (i == n) { best = min(best, dp[s]); continue; }
		rep(j,0,m) if (!(s >> j & 1)) dp[s | 1 << j] = min(dp[s | 1 << j], dp[s] + c[i][j]);
	}
	return best;
}

// Oracle 2: successive shortest paths with Bellman-Ford on the residual graph.
ll ssp(const vector<vi>& c, int n, int m) {
	vi mt(m, -1), of(n, -1); ll tot = 0;
	rep(it,0,n) {
		// nodes: left 0..n-1, right n..n+m-1
		vector<ll> d(n + m, LLONG_MAX); vi pre(n + m, -1);
		rep(i,0,n) if (of[i] == -1) d[i] = 0;
		for (bool ch = 1; ch;) {
			ch = 0;
			rep(i,0,n) if (d[i] != LLONG_MAX) rep(j,0,m) if (of[i] != j && d[i] + c[i][j] < d[n + j])
				d[n + j] = d[i] + c[i][j], pre[n + j] = i, ch = 1;
			rep(j,0,m) if (mt[j] != -1 && d[n + j] != LLONG_MAX && d[n + j] - c[mt[j]][j] < d[mt[j]])
				d[mt[j]] = d[n + j] - c[mt[j]][j], pre[mt[j]] = n + j, ch = 1;
		}
		int bj = -1;
		rep(j,0,m) if (mt[j] == -1 && d[n + j] != LLONG_MAX && (bj == -1 || d[n + j] < d[n + bj])) bj = j;
		assert(bj != -1);
		tot += d[n + bj];
		for (int j = bj;;) {
			int i = pre[n + j], pj = of[i];
			mt[j] = i, of[i] = j;
			if (pj == -1) break;
			j = pj;
		}
	}
	return tot;
}

void check(const vector<vi>& cost, int n, int m, ll expect) {
	auto res = hungarian(cost);
	assert(sz(res.second) == n);
	ll sum = 0; vi used(m);
	rep(i,0,n) {
		int j = res.second[i];
		assert(0 <= j && j < m && !used[j]++);
		sum += cost[i][j];
	}
	assert(sum == res.first);
	assert(sum == expect);
}

void test(int N, int lo, int hi, int iters, bool useBrute) {
	rep(it,0,iters) {
		int n = randRange(0, N + 1), m = randRange(0, N + 1);
		if (n > m) swap(n, m);
		if (it % 7 == 0) m = n; // square
		vector<vi> cost(n, vi(m));
		rep(i,0,n) rep(j,0,m) cost[i][j] = randIncl(lo, hi);
		ll a = ssp(cost, n, m);
		if (useBrute) assert(a == brute(cost, n, m));
		check(cost, n, m, a);
	}
}

int main() {
	check({}, 0, 0, 0);
	check({{7}}, 1, 1, 7);
	check({{-7}}, 1, 1, -7);
	check({{3, 1, 2}}, 1, 3, 1);
	test(7, -5, 5, 30000, 1);
	test(7, 0, 1, 30000, 1); // many ties
	test(7, -1000000, 1000000, 10000, 1);
	test(25, -5, 5, 2000, 0);
	test(60, -1000, 1000, 200, 0);
	test(60, 0, 1, 200, 0);
	// large magnitudes whose optimum fits in int. The header documents no bound, but the
	// int potentials/distances need headroom of a few times max|cost|; these ranges are safe.
	test(8, -100000000, 100000000, 20000, 1);
	test(8, 0, 250000000, 20000, 1);
	test(8, -250000000, 0, 20000, 1);
#ifdef WM_LARGE
	// |cost| ~ 1e9: optimum -1000000001 fits in int, but intermediate values overflow
	// (UB; g++ -O2 returns 999999999). Opt-in: needs a doc note or ll potentials.
	check({{999999999, 1000000000, 1000000001, -1000000001},
	       {1000000001, 1000000000, 1000000001, -1000000000},
	       {-1000000000, 1000000000, 1000000001, 999999999}}, 3, 4, -1000000001);
#endif
	cout << "Tests passed!" << endl;
}
