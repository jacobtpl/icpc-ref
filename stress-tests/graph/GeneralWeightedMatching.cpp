#include "../utilities/template.h"
// from content/contest/template.cpp
#define pb push_back
template<class T> bool ckmax(T &a, const T& b) {return b>a?a=b,1:0;}
template<class T> bool ckmin(T &a, const T& b) {return b<a?a=b,1:0;}

#include "../../content/graph/GeneralWeightedMatching.h"

// Opt-in cases for undocumented preconditions (see the audit findings):
//  -DSZ_EQ_N      N == SZ indexes aux[2*SZ] / floFrom[.][SZ] out of bounds
//                 (visible with -fsanitize=undefined)
//  -DBIG_WEIGHTS  weights >= 2^29 overflow int and give wrong matchings
const int MAXN = 40;
#ifdef SZ_EQ_N
static WeightedMatch<9> W;
#else
static WeightedMatch<MAXN + 1> W;
#endif

mt19937 rng(99);
int ri(int a, int b) { return (int)(rng() % (unsigned)(b - a + 1)) + a; }

// Runs the matcher (the instance is reused on purpose) and validates the
// returned matching against the reported weight and size.
pair<ll,int> run(int n, const vector<vi>& w) {
	W.init(n);
	rep(i,0,n) rep(j,0,i) if (w[i][j]) W.ae(i + 1, j + 1, w[i][j]);
	auto res = W.calc();
	ll ws = 0; int c = 0;
	rep(i,1,n+1) if (int m = W.match[i]) {
		assert(1 <= m && m <= n && m != i && W.match[m] == i);
		assert(w[i-1][m-1] > 0);
		if (m < i) ws += w[i-1][m-1], c++;
	}
	assert(ws == res.first && c == res.second);
	return res;
}

// Max weight matching (cardinality is NOT maximised) by bitmask DP.
ll brute(int n, const vector<vi>& w) {
	vector<ll> dp(1 << n);
	rep(m,1,1<<n) {
		int i = __builtin_ctz(m);
		dp[m] = dp[m ^ (1 << i)];
		rep(j,i+1,n) if ((m >> j & 1) && w[i][j])
			dp[m] = max(dp[m], dp[m ^ (1 << i) ^ (1 << j)] + w[i][j]);
	}
	return dp[(1 << n) - 1];
}

// Max weight bipartite matching (left = [0,a), right = [a,n)), Hungarian.
ll hungarian(int a, int n, const vector<vi>& w) {
	int k = max(a, n - a);
	vector<vector<ll>> c(k + 1, vector<ll>(k + 1));
	rep(i,0,a) rep(j,a,n) c[i + 1][j - a + 1] = -w[i][j];
	vector<ll> u(k + 1), v(k + 1); vi p(k + 1), way(k + 1);
	rep(i,1,k+1) {
		p[0] = i; int j0 = 0;
		vector<ll> minv(k + 1, LLONG_MAX); vector<char> used(k + 1);
		do {
			used[j0] = 1; int i0 = p[j0], j1 = 0; ll delta = LLONG_MAX;
			rep(j,1,k+1) if (!used[j]) {
				ll cur = c[i0][j] - u[i0] - v[j];
				if (cur < minv[j]) minv[j] = cur, way[j] = j0;
				if (minv[j] < delta) delta = minv[j], j1 = j;
			}
			rep(j,0,k+1)
				if (used[j]) u[p[j]] += delta, v[j] -= delta;
				else minv[j] -= delta;
			j0 = j1;
		} while (p[j0]);
		do { int j1 = way[j0]; p[j0] = p[j1]; j0 = j1; } while (j0);
	}
	return v[0];
}

vector<vi> gen(int n, int maxw, bool bip, int a) {
	vector<vi> w(n, vi(n));
	int dens = ri(0, 3) ? ri(0, 100) : 100, mw = ri(0, 1) ? maxw : ri(1, maxw);
	rep(i,0,n) rep(j,0,i) if (ri(0, 99) < dens && (!bip || (j < a) != (i < a)))
		w[i][j] = w[j][i] = ri(1, mw);
	return w;
}

void test(int maxn, int iters, int maxw) {
	rep(it,0,iters) {
		int n = ri(0, maxn);
		auto w = gen(n, maxw, 0, 0);
		assert(run(n, w).first == brute(n, w));
	}
}

void testBipartite(int maxn, int iters, int maxw) {
	rep(it,0,iters) {
		int n = ri(1, maxn), a = ri(0, n);
		auto w = gen(n, maxw, 1, a);
		assert(run(n, w).first == hungarian(a, n, w));
	}
}

int main() {
#ifdef SZ_EQ_N
	test(9, 2000, 5);
#else
	test(1, 10, 5);
	test(3, 20000, 3);
	test(6, 50000, 2);
	test(9, 50000, 5);
	test(9, 20000, 1);
	test(10, 20000, 1000);
	test(12, 3000, (1 << 29) - 1);
	test(16, 200, 10);
	test(18, 30, 1000000);
	testBipartite(MAXN, 2000, 3);
	testBipartite(MAXN, 2000, 500000000);
	{ // complete graph, equal weights; odd cycle; path
		rep(n,1,13) {
			vector<vi> w(n, vi(n, 7)), c(n, vi(n)), p(n, vi(n));
			rep(i,0,n) w[i][i] = 0;
			rep(i,0,n) if (n > 2) c[i][(i+1)%n] = c[(i+1)%n][i] = i + 1;
			rep(i,1,n) p[i][i-1] = p[i-1][i] = 1 + i % 3;
			assert(run(n, w).first == brute(n, w));
			assert(run(n, c).first == brute(n, c));
			assert(run(n, p).first == brute(n, p));
		}
	}
#ifdef BIG_WEIGHTS
	test(9, 20000, (1 << 30) - 1);
#endif
#endif
	cout<<"Tests passed!"<<endl;
}
