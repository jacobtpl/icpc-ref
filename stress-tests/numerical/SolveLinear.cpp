#include "../utilities/template.h"

#include "../../content/numerical/SolveLinear.h"

// Old brute-force test of a copy of the algorithm working mod 3.
namespace modp {

const int mod = 3;
const int nmax = 4, mmax = 4, nmmax = 10;

const int lut[9] = {-4,-2,-3,-1,-100,1,3,2,4};

int modinv(int x) {
	assert(x);
	return x;
	// return lut[x+4];
}

typedef vector<int> vd;

int solveLinear(vector<vd>& A, vd& b, vd& x) {
	int n = sz(A), m = sz(x), rank = 0, br, bc;
	if (n) assert(sz(A[0]) == m);
	vi col(m); iota(all(col), 0);

	rep(i,0,n) {
		int v, bv = -1;
		rep(r,i,n) rep(c,i,m)
			if ((v = A[r][c])) {
				br = r, bc = c, bv = v;
				goto found;
			}
		rep(j,i,n) if (b[j]) return -1;
		break;
found:
		swap(A[i], A[br]);
		swap(b[i], b[br]);
		swap(col[i], col[bc]);
		rep(j,0,n) swap(A[j][i], A[j][bc]);
		bv = modinv(A[i][i]);
		rep(j,i+1,n) {
			int fac = A[j][i] * bv % mod;
			b[j] = (b[j] - fac * b[i]) % mod;
			rep(k,i+1,m) A[j][k] = (A[j][k] - fac*A[i][k]) % mod;
		}
		rank++;
	}

	x.assign(m, 0);
	for (int i = rank; i--;) {
		b[i] = ((b[i] * modinv(A[i][i]) % mod) + mod) % mod;
		x[col[i]] = b[i];
		rep(j,0,i)
			b[j] = (b[j] - A[j][i] * b[i]);
	}
	return rank;
}

template<class F>
void rec(int i, int j, vector<vd>& A, F f) {
	if (i == sz(A)) {
		f();
	}
	else if (j == sz(A[i])) {
		rec(i+1, 0, A, f);
	}
	else {
		rep(v,0,mod) {
			A[i][j] = v;
			rec(i, j+1, A, f);
		}
	}
}

template<class F>
void rec2(int i, vd& A, F f) {
	if (i == sz(A)) f();
	else {
		rep(v,0,mod) {
			A[i] = v;
			rec2(i+1, A, f);
		}
	}
}

void bruteMod() {
	rep(n,0,nmax+1) rep(m,0,mmax+1) {
		int nm = n*m;
		if (nm > nmmax) continue;
		vector<vd> A(n, vd(m));
		vd b(n), x(m), theX(m);
		rec(0, 0, A, [&]() {
			rec2(0, b, [&]() {
				int sols = 0;
				rec2(0, x, [&]() {
					rep(i,0,n) {
						int v = 0;
						rep(j,0,m) v += A[i][j] * x[j];
						if (v % mod != b[i]) return;
					}
					sols++;
					if (sols == 1) theX = x;
				});
				vector<vd> A2 = A;
				vd x2 = x, b2 = b;
				int r = solveLinear(A2, b2, x2);
				if (sols == 0) assert(r == -1);
				else if (sols == 1) assert(r == m);
				else assert(r < m);
				if (sols == 1) assert(x2 == theX);
			});
		});
	}
}
} // namespace modp

// ---- tests of the actual header ----
mt19937 rng(777);
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
	int cnt[3] = {};
	rep(it,0,iters) {
		int n = rnd(0, maxN), m = rnd(0, maxM);
		vector<vector<ll>> A(n, vector<ll>(m)), Ab;
		vector<ll> b(n);
		int mode = rnd(0, 2);
		rep(i,0,n) {
			if (mode && i && rnd(0, 2) == 0) { // linear combination of earlier rows
				int k = rnd(0, i-1), f = rnd(-2, 2);
				rep(j,0,m) A[i][j] = f * A[k][j];
				b[i] = mode == 1 ? f * b[k] : rnd(-C, C);
			} else {
				rep(j,0,m) A[i][j] = rnd(-C, C);
				b[i] = rnd(-C, C);
			}
		}
		Ab = A;
		rep(i,0,n) Ab[i].push_back(b[i]);
		int rk = rankP(A), rk2 = rankP(Ab);
		vector<vd> A2(n, vd(m)); vd b2(n), x(m, 123);
		rep(i,0,n) { rep(j,0,m) A2[i][j] = (double)A[i][j]; b2[i] = (double)b[i]; }
		int r = solveLinear(A2, b2, x);
		if (rk != rk2) { assert(r == -1); cnt[0]++; continue; }
		assert(r == rk); cnt[1 + (rk < m)]++;
		assert(sz(x) == m);
		rep(i,0,n) {
			double s = 0;
			rep(j,0,m) s += (double)A[i][j] * x[j];
			assert(abs(s - (double)b[i]) < 1e-7);
		}
	}
	rep(i,0,3) assert(cnt[i] > iters / 20);
}

void testLarge() {
	// full rank, well conditioned-ish random systems: check residual
	for (int n : {1, 2, 10, 50, 200}) rep(it,0,n <= 50 ? 20 : 2) {
		vector<vd> A(n, vd(n)), A0; vd b(n), b0, x(n);
		rep(i,0,n) { rep(j,0,n) A[i][j] = rnd(-1000000, 1000000) / 1e6; b[i] = rnd(-1000000, 1000000) / 1e6; }
		A0 = A; b0 = b;
		assert(solveLinear(A, b, x) == n);
		rep(i,0,n) {
			double s = 0;
			rep(j,0,n) s += A0[i][j] * x[j];
			assert(abs(s - b0[i]) < 1e-7);
		}
	}
	// underdetermined / overdetermined consistent systems with O(1) entries
	rep(it,0,200) {
		int n = rnd(1, 30), m = rnd(1, 30), r = rnd(1, min(n, m));
		vector<vd> A(n, vd(m)), A0; vd b(n), b0, x(m), x0(m);
		vector<vd> L(n, vd(r)), R(r, vd(m));
		rep(i,0,n) rep(k,0,r) L[i][k] = rnd(-1000, 1000) / 1e3;
		rep(k,0,r) rep(j,0,m) R[k][j] = rnd(-1000, 1000) / 1e3;
		rep(i,0,n) rep(j,0,m) rep(k,0,r) A[i][j] += L[i][k] * R[k][j] / r;
		rep(j,0,m) x0[j] = rnd(-1000, 1000) / 1e3;
		rep(i,0,n) rep(j,0,m) b[i] += A[i][j] * x0[j];
		A0 = A; b0 = b;
		int res = solveLinear(A, b, x);
		assert(res == r);
		rep(i,0,n) {
			double s = 0;
			rep(j,0,m) s += A0[i][j] * x[j];
			assert(abs(s - b0[i]) < 1e-7);
		}
	}
}

void testEdge() {
	vector<vd> A; vd b, x;
	assert(solveLinear(A, b, x) == 0 && x.empty()); // n = m = 0
	x.assign(3, 5);
	assert(solveLinear(A, b, x) == 0 && x == vd(3, 0)); // no equations
	A.assign(2, vd()); b = {0, 0}; x.clear();
	assert(solveLinear(A, b, x) == 0); // no unknowns, consistent
	A.assign(2, vd()); b = {0, 1}; x.clear();
	assert(solveLinear(A, b, x) == -1); // 0 = 1
	A = {{0, 0}, {0, 0}}; b = {0, 0}; x.assign(2, 1);
	assert(solveLinear(A, b, x) == 0 && x == vd(2, 0));
	A = {{2}}; b = {3}; x.assign(1, 0);
	assert(solveLinear(A, b, x) == 1 && abs(x[0] - 1.5) < 1e-12);
#ifdef SOLVELINEAR_ABS_EPS
	// Consistent rank-2 integer system (row 3 = -33 * row 1 + 19 * row 2), entries < 400.
	// eps = 1e-12 is an absolute tolerance, so round-off of ~1e-13 * |A| is treated as "0 = nonzero".
	A = {{25,-26,0},{24,-25,0},{-369,383,0}}; b = {-16,31,1117}; x.assign(3, 0);
	int r = solveLinear(A, b, x);
	cout << "returned " << r << ", expected 2" << endl;
	assert(r == 2);
#endif
}

int main() {
	modp::bruteMod();
	testEdge();
	testSmall(300000, 4, 4, 3);
	testSmall(100000, 6, 6, 2);
	testSmall(100000, 3, 6, 9);
	testLarge();
	cout<<"Tests passed!"<<endl;
}
