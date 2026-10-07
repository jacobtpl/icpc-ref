#include "../utilities/template.h"

#include "../../content/numerical/SolveLinear.h"

// SolveLinear2.h is a diff against SolveLinear.h rather than a compilable header, so this is
// SolveLinear.h's solveLinear with exactly the changes described in SolveLinear2.h applied.
const double undefined = 1e300;
int solveLinear2(vector<vd>& A, vd& b, vd& x) {
	int n = sz(A), m = sz(x), rank = 0, br, bc;
	if (n) assert(sz(A[0]) == m);
	vi col(m); iota(all(col), 0);

	rep(i,0,n) {
		double v, bv = 0;
		rep(r,i,n) rep(c,i,m)
			if ((v = fabs(A[r][c])) > bv)
				br = r, bc = c, bv = v;
		if (bv <= eps) {
			rep(j,i,n) if (fabs(b[j]) > eps) return -1;
			break;
		}
		swap(A[i], A[br]);
		swap(b[i], b[br]);
		swap(col[i], col[bc]);
		rep(j,0,n) swap(A[j][i], A[j][bc]);
		bv = 1/A[i][i];
		rep(j,0,n) if (j != i) { // instead of rep(j,i+1,n)
			double fac = A[j][i] * bv;
			b[j] -= fac * b[i];
			rep(k,i+1,m) A[j][k] -= fac*A[i][k];
		}
		rank++;
	}

	x.assign(m, undefined);
	rep(i,0,rank) {
		rep(j,rank,m) if (fabs(A[i][j]) > eps) goto fail;
		x[col[i]] = b[i] / A[i][i];
	fail:; }
	return rank;
}

// Old timing-only test of a copy.

enum { YES, NO, MULT };
int solve_linear(vector<vd>& A, vd& b, vd& x) {
	int n = sz(A), m = sz(x), br = -1, bc = -1;
	vi col(m); iota(all(col), 0);

	rep(i,0,n) {
		double v, bv = -1;
		rep(r,i,n) rep(c,i,m)
			if ((v = fabs(A[r][c])) > bv)
				br = r, bc = c, bv = v;
		if (bv <= eps) {
			rep(j,i,n) if (fabs(b[j]) > eps) return NO;
			if (i == m) break;
			return MULT;
		}
		swap(A[i], A[br]);
		swap(b[i], b[br]);
		swap(col[i], col[bc]);
		rep(j,0,n) swap(A[j][i], A[j][bc]);
		bv = 1/A[i][i];
		rep(j,i+1,n) {
			double fac = A[j][i] * bv;
			b[j] -= fac * b[i];
			rep(k,i+1,m) A[j][k] -= fac*A[i][k];
		}
	}
	if (n < m) return MULT;

	for (int i = m; i--;) {
		x[col[i]] = (b[i] /= A[i][i]);
		rep(j,0,i)
			b[j] -= A[j][i] * b[i];
	}
	return YES;
}

void oldTest() {
	const int n = 1000;
	vector<vd> A(n, vd(n));
	rep(i,0,n) rep(j,0,n) A[i][j] = rand() * 1000.0 / RAND_MAX;
	vd x(n), b(n);
	rep(i,0,n) b[i] = rand() * 1000.0 / RAND_MAX;
	int r = solve_linear(A, b, x);
	assert(r == 0);
	// cout << r << endl;
	// rep(i,0,n) cout << x[i] << ' ';
	// cout << endl;
}

mt19937 rng(4242);
int rnd(int lo, int hi) { return uniform_int_distribution<int>(lo, hi)(rng); }

const ll P = 1000000007;
ll mpow(ll a, ll e) { ll r = 1; for (a %= P; e; e /= 2, a = a * a % P) if (e & 1) r = r * a % P; return r; }
// rank mod P; exact for small integer matrices (all minors are far below P)
int rankP(vector<vector<ll>> a) {
	int n = sz(a), m = n ? sz(a[0]) : 0, r = 0;
	for (auto& row : a) for (auto& v : row) v = (v % P + P) % P;
	rep(c,0,m) {
		int p = r;
		while (p < n && !a[p][c]) p++;
		if (p == n) continue;
		swap(a[p], a[r]);
		ll iv = mpow(a[r][c], P - 2);
		rep(i,r+1,n) {
			ll f = a[i][c] * iv % P;
			rep(k,c,m) a[i][k] = (a[i][k] - f * a[r][k] % P + P) % P;
		}
		if (++r == n) break;
	}
	return r;
}

void testSmall(int iters, int maxN, int maxM, int C) {
	int cnt[4] = {};
	rep(it,0,iters) {
		int n = rnd(0, maxN), m = rnd(0, maxM);
		vector<vector<ll>> A(n, vector<ll>(m)), Ab;
		vector<ll> b(n);
		int mode = rnd(0, 2);
		rep(i,0,n) {
			if (mode && i && rnd(0, 2) == 0) {
				int k = rnd(0, i-1), f = rnd(-2, 2);
				rep(j,0,m) A[i][j] = f * A[k][j];
				b[i] = mode == 1 ? f * b[k] : rnd(-C, C);
			} else {
				rep(j,0,m) A[i][j] = rnd(0, 2) ? rnd(-C, C) : 0;
				b[i] = rnd(-C, C);
			}
		}
		Ab = A;
		rep(i,0,n) Ab[i].push_back(b[i]);
		int rk = rankP(A), rk2 = rankP(Ab);
		vector<vd> A2(n, vd(m)), A3; vd b2(n), b3, x(m, 123), y(m);
		rep(i,0,n) { rep(j,0,m) A2[i][j] = (double)A[i][j]; b2[i] = (double)b[i]; }
		A3 = A2; b3 = b2;
		int r = solveLinear2(A2, b2, x);
		if (rk != rk2) { assert(r == -1); cnt[0]++; continue; }
		assert(r == rk);
		assert(solveLinear(A3, b3, y) == rk); // some solution
		rep(j,0,m) {
			// x_j is uniquely determined iff e_j lies in the row space of A
			auto Ae = A; Ae.emplace_back(m); Ae.back()[j] = 1;
			bool uniq = rankP(Ae) == rk;
			cnt[1 + uniq]++;
			if (uniq) assert(x[j] != undefined && abs(x[j] - y[j]) < 1e-7);
			else assert(x[j] == undefined);
		}
		cnt[3] += rk < m && count(all(x), undefined) < m; // mixed case
	}
	rep(i,0,4) assert(cnt[i] > iters / 20);
}

int main() {
	oldTest();
	{
		vector<vd> A; vd b, x;
		assert(solveLinear2(A, b, x) == 0 && x.empty());
		x.assign(2, 0);
		assert(solveLinear2(A, b, x) == 0 && x == vd(2, undefined));
		// x0 + x1 = 2, x2 = 3: only x2 is determined
		A = {{1,1,0},{0,0,1}}; b = {2,3}; x.assign(3, 0);
		assert(solveLinear2(A, b, x) == 2 && x[0] == undefined && x[1] == undefined && abs(x[2] - 3) < 1e-12);
	}
	testSmall(300000, 4, 4, 3);
	testSmall(100000, 6, 6, 2);
	testSmall(100000, 3, 6, 9);
	cout<<"Tests passed!"<<endl;
}
