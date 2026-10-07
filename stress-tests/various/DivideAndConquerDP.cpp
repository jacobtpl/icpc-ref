#include "../utilities/template.h"

// The header is a template to be modified ("Modify at will"); out of the box it needs `dp`
// (cost table, f(ind, k) = dp[ind][k]) and `res` (output) to be defined by the user, and
// uses lo(ind) = 0, hi(ind) = ind, i.e. a[i] = min_{0 <= k < i} dp[i][k].
// Compile with -DBENCH for timings.
ll evals;
int mode; // 0: explicit matrix, 1: implicit Monge cost, 2: constant cost
vector<vector<ll>> mat;
vector<ll> pre, g;
struct Row {
	int i;
	ll operator[](int k) const {
		evals++;
		if (mode == 0) return mat[i][k];
		if (mode == 2) return 5;
		ll d = pre[i] - pre[k];
		return g[k] + d * d;
	}
};
struct { Row operator[](int i) const { return {i}; } } dp;
vector<pii> res;

#include "../../content/various/DivideAndConquerDP.h"

mt19937 rng(2024);
int rnd() { return (int)(rng() >> 1); }

// minimal optimal k and the optimum for row i, by brute force
pii brute(int i) {
	pair<ll, int> best(LLONG_MAX, -1);
	rep(k,0,i) best = min(best, make_pair(mat[i][k], k));
	return pii(best.second, (int)best.first);
}

// Checks solve(L, R) on the current matrix (which must have monotone minimal argmins).
void check(int n) {
	vector<pii> want(n);
	rep(i,1,n) want[i] = brute(i);
	rep(i,2,n) assert(want[i-1].first <= want[i].first);
	rep(L,0,min(n,3)+1) for (int R : {n, n - 1, L + 1, L}) {
		if (R < L || R > n) continue;
		res.assign(n, pii(-7, -7));
		DP().solve(L, R);
		rep(i,0,n) {
			if (i < L || i >= R) assert(res[i] == pii(-7, -7)); // untouched
			else if (i > 0) assert(res[i] == want[i]);
			else assert(res[i].first == INT_MIN); // no candidates for i = 0
		}
	}
}

int main() {
	mode = 0;
	// Exhaustive: all cost tables with entries in {0,1,2} for n <= 6 whose minimal argmin is monotone.
	rep(n,1,7) {
		int cells = n * (n - 1) / 2, total = 1;
		rep(i,0,cells) total *= 3;
		mat.assign(n, vector<ll>(n));
		rep(code,0,total) {
			int c = code, prevOpt = 0; bool ok = 1;
			rep(i,1,n) {
				rep(k,0,i) mat[i][k] = c % 3, c /= 3;
				int o = (int)(min_element(mat[i].begin(), mat[i].begin() + i) - mat[i].begin());
				if (o < prevOpt) { ok = 0; break; }
				prevOpt = o;
			}
			if (!ok) continue;
			if (n == 6) { // only solve(1, n) for the largest size, for speed
				res.assign(n, pii());
				DP().solve(1, n);
				rep(i,1,n) assert(res[i] == brute(i));
			} else check(n);
		}
	}
	// Random tables with a planted monotone minimal argmin, many ties, negative values.
	rep(it,0,100000) {
		int n = rnd() % 25 + 1, range = rnd() % 3 ? 3 : 1000000;
		mat.assign(n, vector<ll>(n));
		int o = 0;
		rep(i,1,n) {
			o = min(i - 1, o + (rnd() % 3 == 0 ? rnd() % 4 : 0));
			int mn = rnd() % (2 * range) - range;
			rep(k,0,i) mat[i][k] = mn + (k < o ? 1 : 0) + (k == o ? 0 : rnd() % range);
		}
		check(n);
	}
	// Monge costs: a[i] = min_k g[k] + (pre[i] - pre[k])^2; compare against the O(n^2) brute force.
	rep(it,0,20000) {
		int n = rnd() % 60 + 1;
		pre.assign(n, 0); g.assign(n, 0);
		rep(i,1,n) pre[i] = pre[i-1] + rnd() % 20;
		rep(i,0,n) g[i] = rnd() % 2000 - 1000;
		mode = 1;
		mat.assign(n, vector<ll>(n));
		rep(i,0,n) rep(k,0,i) mat[i][k] = dp[i][k];
		res.assign(n, pii());
		evals = 0;
		DP().solve(1, n);
		rep(i,1,n) assert(res[i] == brute(i));
		int lg = 1; while ((1 << lg) < n) lg++;
		assert(evals <= (ll)(2 * n + 2) * (lg + 1)); // O(N log N) evaluations
		mode = 0;
		check(n);
	}

#ifdef BENCH
	for (int md : {1, 2}) for (int n : {100000, 1000000, 2000000}) {
		mode = md;
		pre.assign(n, 0); g.assign(n, 0);
		rep(i,1,n) pre[i] = pre[i-1] + rnd() % 20;
		rep(i,0,n) g[i] = rnd() % 2000000 - 1000000;
		res.assign(n, pii());
		evals = 0;
		auto t0 = chrono::steady_clock::now();
		DP().solve(1, n);
		auto t1 = chrono::steady_clock::now();
		cerr << (md == 1 ? "monge" : "constant") << " cost N=" << n << ": " << evals << " evals of f, "
			<< chrono::duration<double>(t1 - t0).count() << " s" << endl;
	}
#endif
	cout<<"Tests passed!"<<endl;
}
