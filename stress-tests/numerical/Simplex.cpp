#include "../utilities/template.h"

#include "../../content/numerical/Simplex.h"

mt19937 rng(2024);
int rnd(int lo, int hi) { return uniform_int_distribution<int>(lo, hi)(rng); }

// ---- exact brute-force oracle: enumerate vertices with integer arithmetic ----
ll det(vector<vector<ll>> a) {
	int n = sz(a);
	if (n == 0) return 1;
	if (n == 1) return a[0][0];
	ll r = 0;
	rep(c,0,n) {
		vector<vector<ll>> b;
		rep(i,1,n) { b.emplace_back(); rep(j,0,n) if (j != c) b.back().push_back(a[i][j]); }
		r += (c % 2 ? -1 : 1) * a[0][c] * det(b);
	}
	return r;
}

// maximise c.x over {G x <= h} (G already contains x >= 0 rows); the region is pointed.
// returns (feasible, value as num/den with den > 0)
struct Res { bool feas; ll num, den; };
Res vertexOpt(const vector<vector<ll>>& G, const vector<ll>& h, const vector<ll>& c) {
	int k = sz(G), n = sz(c);
	Res best{false, 0, 1};
	vi pick(k, 0);
	fill(pick.begin(), pick.begin() + n, 1);
	do {
		vector<vector<ll>> M; vector<ll> rhs;
		rep(i,0,k) if (pick[i]) M.push_back(G[i]), rhs.push_back(h[i]);
		ll d = det(M);
		if (d == 0) continue;
		vector<ll> xn(n); // x[j] = xn[j] / d
		rep(j,0,n) {
			auto M2 = M;
			rep(i,0,n) M2[i][j] = rhs[i];
			xn[j] = det(M2);
		}
		if (d < 0) { d = -d; for (auto& v : xn) v = -v; }
		bool ok = true;
		rep(i,0,k) {
			ll s = 0;
			rep(j,0,n) s += G[i][j] * xn[j];
			if (s > h[i] * d) { ok = false; break; }
		}
		if (!ok) continue;
		ll num = 0;
		rep(j,0,n) num += c[j] * xn[j];
		if (!best.feas || num * best.den > best.num * d) best = {true, num, d};
	} while (prev_permutation(all(pick)));
	return best;
}

void checkFeasible(const vvd& A, const vd& b, const vd& x, double tol) {
	int m = sz(b), n = sz(x);
	rep(j,0,n) assert(x[j] >= -tol);
	rep(i,0,m) {
		double s = 0;
		rep(j,0,n) s += A[i][j] * x[j];
		assert(s <= b[i] + tol);
	}
}

void testExact(int iters, int maxN, int maxM, int C) {
	int cntInf = 0, cntUnb = 0, cntOpt = 0;
	rep(it,0,iters) {
		int n = rnd(1, maxN), m = rnd(0, maxM);
		vector<vector<ll>> G; vector<ll> h, c(n);
		vvd A(m, vd(n)); vd b(m), cc(n), x;
		int zeroRhs = rnd(0, 3) == 0; // degenerate vertex at the origin
		rep(i,0,m) {
			G.emplace_back(n);
			rep(j,0,n) A[i][j] = (double)(G[i][j] = rnd(-C, C));
			h.push_back(zeroRhs ? 0 : rnd(-C, C)); b[i] = (double)h[i];
		}
		rep(j,0,n) {
			cc[j] = (double)(c[j] = rnd(-C, C));
			G.emplace_back(n); G.back()[j] = -1; h.push_back(0);
		}
		Res r = vertexOpt(G, h, c);
		bool unb = false;
		if (r.feas) {
			// recession cone: G d <= 0, sum d <= 1
			auto G2 = G; vector<ll> h2(sz(G), 0);
			G2.emplace_back(n, 1); h2.push_back(1);
			Res ray = vertexOpt(G2, h2, c);
			assert(ray.feas);
			unb = ray.num > 0;
		}
		T v = LPSolver(A, b, cc).solve(x);
		if (!r.feas) { assert(v == -inf); cntInf++; }
		else if (unb) {
			assert(v == inf); cntUnb++;
			assert(sz(x) == n); checkFeasible(A, b, x, 1e-6);
		} else {
			assert(v != inf && v != -inf); cntOpt++;
			assert(abs(v - (double)r.num / (double)r.den) < 1e-6);
			assert(sz(x) == n); checkFeasible(A, b, x, 1e-6);
			double s = 0;
			rep(j,0,n) s += cc[j] * x[j];
			assert(abs(s - v) < 1e-6);
		}
	}
	assert(cntInf > iters / 50 && cntUnb > iters / 50 && cntOpt > iters / 50);
}

// larger instances: certify optimality through the dual (weak duality).
void testDual(int iters, int maxN, int maxM, int C, bool origin) {
	int both = 0;
	rep(it,0,iters) {
		int n = rnd(1, maxN), m = rnd(1, maxM);
		vvd A(m, vd(n)), At(n, vd(m)); vd b(m), c(n), x, y, nb(m), nc(n);
		rep(i,0,m) rep(j,0,n) {
			A[i][j] = rnd(0, 2) ? rnd(-C, C) : 0;
			At[j][i] = -A[i][j];
		}
		// bias towards feasible + bounded: positive rows keep it bounded
		rep(j,0,n) A[rnd(0, m-1)][j] = rnd(1, C), c[j] = rnd(-C, C);
		rep(i,0,m) rep(j,0,n) At[j][i] = -A[i][j];
		rep(i,0,m) b[i] = origin ? rnd(0, C) : rnd(-C / 4, 3 * C);
		rep(i,0,m) nb[i] = -b[i];
		rep(j,0,n) nc[j] = -c[j];
		T p = LPSolver(A, b, c).solve(x);
		T d = -LPSolver(At, nc, nb).solve(y); // min b.y, A^T y >= c, y >= 0
		if (p == -inf) { assert(d == -inf || d == inf); continue; }
		checkFeasible(A, b, x, 1e-6);
		if (p == inf) { assert(d == inf); continue; } // dual infeasible: -(-inf)
		both++;
		assert(abs(p - d) < 1e-6 * max(1.0, abs(p)));
		checkFeasible(At, nc, y, 1e-6);
		double s = 0;
		rep(j,0,n) s += c[j] * x[j];
		assert(abs(s - p) < 1e-6 * max(1.0, abs(p)));
	}
	assert(both > iters / 10);
}

void testEdge() {
	vd x;
	{ // usage example from the header
		vvd A = {{1,-1}, {-1,1}, {-1,-2}};
		vd b = {1,1,-4}, c = {-1,-1};
		T val = LPSolver(A, b, c).solve(x);
		assert(abs(val + 7.0 / 3) < 1e-9);
		assert(abs(x[0] - 2.0 / 3) < 1e-9 && abs(x[1] - 5.0 / 3) < 1e-9);
	}
	{ // no constraints
		assert(LPSolver({}, {}, {-1, 0}).solve(x) == 0 && x == vd({0, 0}));
		assert(LPSolver({}, {}, {-1, 1}).solve(x) == inf && sz(x) == 2);
	}
	{ // x <= -1 infeasible; x <= 0 forces x = 0
		assert(LPSolver({{1}}, {-1}, {1}).solve(x) == -inf);
		assert(LPSolver({{1}}, {0}, {1}).solve(x) == 0);
		assert(LPSolver({{-1}}, {-5}, {-1}).solve(x) == -5 && abs(x[0] - 5) < 1e-9);
	}
	{ // duplicated / redundant constraints, equality via two inequalities
		vvd A = {{1,1},{1,1},{-1,-1},{1,0},{1,0}};
		vd b = {4,4,-4,3,3}, c = {2,1};
		assert(abs(LPSolver(A, b, c).solve(x) - 7) < 1e-9);
		assert(abs(x[0] - 3) < 1e-9 && abs(x[1] - 1) < 1e-9);
	}
	{ // Klee-Minty cube: optimum 5^n (exponential for Dantzig's rule, must still be right)
		for (int n : {1, 2, 5, 10, 13}) {
			vvd A(n, vd(n)); vd b(n), c(n);
			rep(i,0,n) {
				rep(j,0,i) A[i][j] = pow(2.0, i - j + 1);
				A[i][i] = 1; b[i] = pow(5.0, i + 1); c[i] = pow(2.0, n - 1 - i);
			}
			T v = LPSolver(A, b, c).solve(x);
			assert(abs(v - pow(5.0, n)) < 1e-6 * pow(5.0, n));
		}
	}
}

int main() {
	testEdge();
	testExact(60000, 3, 4, 3);
	testExact(20000, 2, 6, 5);
	testExact(4000, 4, 4, 2);
	testExact(20000, 3, 5, 1);
	testDual(3000, 8, 8, 10, false);
	testDual(3000, 8, 8, 10, true);
	testDual(300, 40, 40, 20, false);
	testDual(300, 40, 40, 20, true);
	cout<<"Tests passed!"<<endl;
}
