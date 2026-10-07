#include "../utilities/template.h"

// From content/contest/template.cpp, which Gaussian.h relies on.
bool ckmax(auto &a, auto const& b) {return b>a?a=b,1:0;}

// Gaussian.h declares a getRow overload for `mint` even when T = double,
// so some `mint` with comparison against 0 has to exist. The notebook's own
// mint (number-theory/SiyongModular.h) has no operator== / operator!=, so
// the header does not compile with it: build with -DGAUSS_SIYONG_MINT to see.
#ifdef GAUSS_SIYONG_MINT
#include "../../content/number-theory/SiyongModular.h"
#else
struct mint {
	int v;
	bool operator!=(int o) const { return v != o; }
};
#endif

#include "../../content/numerical/Gaussian.h"

mt19937_64 rng(8675309);
int rnd(int lo, int hi) { return lo + (int)(rng() % (unsigned)(hi - lo + 1)); }

// Exact rationals for the oracle.
struct Fr {
	ll p, q;
	Fr(ll a = 0, ll b = 1) {
		ll g = __gcd(abs(a), abs(b));
		if (g == 0) g = 1;
		if (b < 0) g = -g;
		p = a / g; q = b / g;
	}
	Fr operator*(Fr o) const { return Fr(p * o.p, q * o.q); }
	Fr operator-(Fr o) const { return Fr(p * o.q - o.p * q, q * o.q); }
	Fr operator/(Fr o) const { return Fr(p * o.q, q * o.p); }
	bool zero() const { return p == 0; }
	double val() const { return (double)p / (double)q; }
};

// Reduced row echelon form (unique), rank and determinant, exactly.
struct Res { vector<vector<Fr>> m; int rank; Fr det; };
Res exact(const vector<vi>& a) {
	int R = sz(a), C = R ? sz(a[0]) : 0, r = 0;
	Res res;
	res.det = Fr(1);
	res.m.assign(R, vector<Fr>(C));
	auto& m = res.m;
	rep(i,0,R) rep(j,0,C) m[i][j] = Fr(a[i][j]);
	rep(c,0,C) {
		int p = -1;
		rep(i,r,R) if (!m[i][c].zero()) { p = i; break; }
		if (p == -1) { res.det = Fr(0); continue; }
		if (p != r) swap(m[p], m[r]), res.det = res.det * Fr(-1);
		Fr piv = m[r][c];
		res.det = res.det * piv;
		rep(k,0,C) m[r][k] = m[r][k] / piv;
		rep(i,0,R) if (i != r && !m[i][c].zero()) {
			Fr v = m[i][c];
			rep(k,0,C) m[i][k] = m[i][k] - v * m[r][k];
		}
		r++;
	}
	res.rank = r;
	if (R != C) res.det = Fr(0);
	return res;
}

void check(const vector<vi>& a) {
	int R = sz(a), C = R ? sz(a[0]) : 0;
	vector<vector<double>> m(R);
	rep(i,0,R) m[i].assign(all(a[i]));
	auto [prod, rank] = gauss(m);
	Res want = exact(a);
	assert(rank == want.rank);
	assert(sz(m) == R);
	rep(i,0,R) {
		assert(sz(m[i]) == C);
		rep(j,0,C) assert(fabs(m[i][j] - want.m[i][j].val()) < 1e-7);
	}
	if (R == C) { // the first component is the determinant
		double d = want.det.val();
		assert(fabs(prod - d) <= 1e-9 * max(1.0, fabs(d)));
		if (rank < R) assert(prod == 0);
	}
}

int main() {
	{ // empty matrix, and rows without columns
		vector<vector<double>> e;
		auto r = gauss(e);
		assert(r.first == 1 && r.second == 0);
		vector<vector<double>> e2(3);
		r = gauss(e2);
		assert(r.second == 0);
	}
	check({{0}}); check({{5}}); check({{-2}});
	check({{0, 0}, {0, 0}});
	check({{0, 1}, {1, 0}});
	check({{1, 2}, {2, 4}});
	check({{0, 0, 1}, {0, 2, 0}, {3, 0, 0}});
	rep(it,0,300000) {
		int R = rnd(1, 6), C = rnd(1, 6), lim = rnd(1, 5), type = rnd(0, 3);
		if (rnd(0, 1)) C = R;
		vector<vi> a(R, vi(C));
		if (type <= 1) {
			for (auto& r : a) for (int& x : r) x = rnd(-lim, lim);
			if (type == 1) for (auto& r : a) for (int& x : r) if (rnd(0, 2)) x = 0;
		} else { // rank at most k: product of R x k and k x C
			int k = rnd(0, min(R, C)), l = type == 2 ? 2 : 1;
			vector<vi> u(R, vi(k)), v(k, vi(C));
			for (auto& r : u) for (int& x : r) x = rnd(-l, l);
			for (auto& r : v) for (int& x : r) x = rnd(-l, l);
			rep(i,0,R) rep(j,0,C) rep(t,0,k) a[i][j] += u[i][t] * v[t][j];
		}
		check(a);
	}
	// Larger well-conditioned system: [A | A x] reduces to [I | x].
	rep(it,0,20) {
		int n = rnd(50, 120);
		vector<vector<double>> m(n, vector<double>(n + 1));
		vector<double> x(n);
		for (double& v : x) v = rnd(-1000, 1000) / 100.0;
		rep(i,0,n) {
			rep(j,0,n) m[i][j] = rnd(-1000, 1000) / 1000.0;
			m[i][i] += n; // diagonally dominant
			rep(j,0,n) m[i][n] += m[i][j] * x[j];
		}
		auto r = gauss(m);
		assert(r.second == n);
		rep(i,0,n) {
			rep(j,0,n) assert(fabs(m[i][j] - (i == j)) < 1e-9);
			assert(fabs(m[i][n] - x[i]) < 1e-7);
		}
	}
	cout<<"Tests passed!"<<endl;
}
