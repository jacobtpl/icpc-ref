#include "../utilities/template.h"
#include "../../content/number-theory/ModInt.h"
#include "../../content/numerical/Determinant.h"
#include "../../content/numerical/Gaussian.h"
#include "../../content/numerical/MatrixInverse.h"

template<class M> M bruteDet(const vector<vector<M>>& a) {
	int n = sz(a); vi perm(n); iota(all(perm), 0);
	M out = 0;
	do {
		M term = 1;
		rep(i,0,n) {
			term *= a[i][perm[i]];
			rep(j,0,i) if (perm[j] > perm[i]) term = -term;
		}
		out += term;
	} while (next_permutation(all(perm)));
	return out;
}

template<int P> void test(int iters, mt19937& rng) {
	using M = ModInt<P>;
	using Mat = vector<vector<M>>;
	Mat empty;
	assert(det(empty) == 1 && matInv(empty) == 0);
	rep(it,0,iters) {
		int n = (int)(rng() % 6) + 1;
		Mat a(n, vector<M>(n));
		// Small values and zeros give many singular matrices and zero pivots.
		for (auto& row : a) for (auto& v : row)
			v = rng() % 3 ? rng() % min(P, 4) : rng() % P;
		if (rng() % 4 == 0 && n > 1) a[n-1] = a[0];
		M want = bruteDet(a);
		Mat c = a;
		assert(det(c) == want);
		c = a; auto [d, rank] = gauss(c);
		assert(d == want && (rank == n) == (want != 0));
		c = a; int r = matInv(c);
		assert((r == n) == (want != 0));
		if (r == n) rep(i,0,n) rep(j,0,n) {
			M s = 0;
			rep(k,0,n) s += a[i][k] * c[k][j];
			assert(s == M(i == j));
		}
	}
}

int main() {
	mt19937 rng(2024);
	test<2>(20000, rng); test<3>(20000, rng); test<7>(20000, rng);
	test<998244353>(20000, rng);
	cout << "Tests passed!" << endl;
}
