#include "../utilities/template.h"

// MatrixInverse.h says "using T = double; // or mint" and ships the
// mint versions of ABS/ISZERO as comments. The header cannot be
// switched from the outside, so this test instantiates the very same
// code with T = mint by renaming `double` while including it; fabs on
// a mint returns its value, which makes the default macros equivalent
// to the commented-out ones (ABS(x) = x.v, ISZERO(x) = (x.v == 0)).
#include "../../content/number-theory/SiyongModular.h"
int fabs(mint x) { return x.v; }
#define double mint
#include "../../content/numerical/MatrixInverse.h"
#undef double

typedef vector<vector<mint>> vvm;
mt19937_64 rng(8088);
int rnd(int lo, int hi) { return uniform_int_distribution<int>(lo, hi)(rng); }

int rankOracle(vvm a) {
	int n = sz(a), r = 0;
	rep(c,0,n) {
		int p = r;
		while (p < n && !a[p][c].v) p++;
		if (p == n) continue;
		swap(a[r], a[p]);
		rep(j,r+1,n) {
			mint f = a[j][c] / a[r][c];
			rep(k,c,n) a[j][k] -= f * a[r][k];
		}
		r++;
	}
	return r;
}

void check(const vvm& m) {
	int n = sz(m);
	vvm inv = m;
	int rk = matInv(inv);
	assert(rk == rankOracle(m));
	if (rk < n) return;
	rep(i,0,n) rep(j,0,n) {
		mint s1, s2;
		rep(k,0,n) s1 += m[i][k] * inv[k][j], s2 += inv[i][k] * m[k][j];
		assert(s1.v == (i == j) && s2.v == (i == j));
	}
}

vvm randMat(int n, int lo, int hi) {
	vvm m(n, vector<mint>(n));
	for (auto& r : m) for (auto& x : r) x = mint(rnd(lo, hi));
	return m;
}

int main() {
	check(vvm());
	// exhaustive 0/1/2 matrices for n <= 3
	rep(n,1,4) {
		int total = 1;
		rep(i,0,n*n) total *= 3;
		rep(mask,0,total) {
			vvm m(n, vector<mint>(n));
			int x = mask;
			rep(i,0,n) rep(j,0,n) m[i][j] = mint(x % 3), x /= 3;
			check(m);
		}
	}
	rep(it,0,30000) {
		int n = rnd(1, it % 100 ? 6 : 30);
		vvm m;
		switch (it % 5) {
		case 0: m = randMat(n, 0, MOD - 1); break;
		case 1: m = randMat(n, 0, 1); break;
		case 2: m = randMat(n, MOD - 3, MOD - 1); break;
		case 3:
			m = randMat(n, 0, MOD - 1);
			for (auto& r : m) for (auto& x : r) if (rnd(0, 2)) x = mint(0);
			break;
		default: {
			int r = rnd(0, n);
			vvm A = randMat(n, 0, MOD - 1), B = randMat(n, 0, MOD - 1);
			m.assign(n, vector<mint>(n));
			rep(i,0,n) rep(j,0,n) rep(k,0,r) m[i][j] += A[i][k] * B[k][j];
		}
		}
		check(m);
	}
	check(randMat(120, 0, MOD - 1));
	cout<<"Tests passed!"<<endl;
}
