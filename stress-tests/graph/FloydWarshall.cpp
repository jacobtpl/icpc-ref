#include "../utilities/template.h"

#include "../../content/graph/FloydWarshall.h"

mt19937_64 rng(4242);
ll rl(ll a, ll b) { return (ll)(rng() % (unsigned long long)(b - a + 1)) + a; }

// Oracle: Bellman-Ford from every source on the edge list, then everything
// reachable from a still-relaxable vertex is -inf.
vector<vector<ll>> oracle(const vector<vector<ll>>& m) {
	int n = sz(m);
	vector<vector<ll>> res(n, vector<ll>(n, inf));
	rep(s,0,n) {
		vector<__int128> d(n); vi reach(n), neg(n);
		reach[s] = 1; d[s] = 0;
		rep(round,0,n) rep(i,0,n) if (reach[i]) rep(j,0,n) if (m[i][j] != inf)
			if (!reach[j] || d[i] + m[i][j] < d[j])
				reach[j] = 1, d[j] = d[i] + m[i][j];
		rep(i,0,n) if (reach[i]) rep(j,0,n) if (m[i][j] != inf)
			if (d[i] + m[i][j] < d[j]) neg[j] = 1;
		rep(round,0,n) rep(i,0,n) if (neg[i]) rep(j,0,n)
			if (m[i][j] != inf) neg[j] = 1;
		rep(i,0,n) if (reach[i]) res[s][i] = neg[i] ? -inf : (ll)d[i];
	}
	return res;
}

void test(int maxn, int iters, ll lo, ll hi, bool loops) {
	rep(it,0,iters) {
		int n = (int)rl(0, maxn), dens = (int)rl(0, 100);
		ll l = rl(0, 3) ? lo : 0; // sometimes non-negative only
		vector<vector<ll>> m(n, vector<ll>(n, inf));
		rep(i,0,n) rep(j,0,n) if ((i != j || loops) && rl(0, 99) < dens)
			m[i][j] = rl(l, hi);
		auto want = oracle(m);
		floydWarshall(m);
		assert(m == want);
	}
}

int main() {
	test(1, 100, -5, 5, 1);
	test(4, 100000, -3, 6, 1);
	test(7, 50000, -2, 10, 0);
	test(8, 30000, -10, 10, 1);
	test(8, 20000, 0, 3, 0);
	test(25, 500, -1, 30, 0);
	test(25, 500, -(ll)1e15, (ll)1e15, 1); // big weights, n*w fits in ll
	test(60, 30, -(ll)1e16, (ll)1e16, 1); // negative cycles blow up to -inf
	test(60, 30, -1, 100, 0);
	{ // long negative cycle: distances would underflow without clamping
		int n = 120;
		vector<vector<ll>> m(n, vector<ll>(n, inf));
		rep(i,0,n) m[i][(i+1)%n] = -(ll)1e17;
		auto want = oracle(m);
		floydWarshall(m);
		assert(m == want);
	}
	cout<<"Tests passed!"<<endl;
}
