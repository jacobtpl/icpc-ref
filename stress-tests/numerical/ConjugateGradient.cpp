#include "../utilities/template.h"

#include "../../content/numerical/ConjugateGradient.h"
#include "../../content/numerical/SolveLinear.h"

template<class T> void testDense(int n) {
	// A = M^T M + I, random b
	vector<vd> M(n, vd(n)), A(n, vd(n));
	rep(i,0,n) rep(j,0,n) M[i][j] = rand() % 21 - 10;
	rep(i,0,n) rep(j,0,n) {
		rep(k,0,n) A[i][j] += M[k][i] * M[k][j];
		if (i == j) A[i][j] += 1;
	}
	vector<T> b(n), dg(n);
	rep(i,0,n) b[i] = rand() % 201 - 100, dg[i] = (T)A[i][i];
	bool ones = rand() % 2;
	auto x = conjGrad([&](const vector<T>& v, vector<T>& y) {
		rep(i,0,n) { y[i] = 0; rep(j,0,n) y[i] += (T)A[i][j] * v[j]; }
	}, b, ones ? vector<T>(n, 1) : dg);
	vd bd(all(b)), xd(n);
	auto A2 = A;
	assert(solveLinear(A2, bd, xd) == n);
	rep(i,0,n) assert(abs((double)x[i] - xd[i]) <= 1e-6 * max(1.0, abs(xd[i])));
}

// connected weighted multigraph, node 0 grounded; answer = sum b[v] * x[v]
template<class T> void testLaplacian(int n, int m) {
	vi U, V, D;
	rep(i,1,n) U.push_back(rand() % i), V.push_back(i), D.push_back(rand() % 100 + 1);
	rep(i,0,m) {
		int u = rand() % n, v = rand() % n;
		if (u != v) U.push_back(u), V.push_back(v), D.push_back(rand() % 100 + 1);
	}
	vector<T> b(n), dg(n);
	vector<vd> A(n - 1, vd(n - 1));
	rep(i,1,n) b[i] = rand() % 101;
	rep(e,0,sz(U)) {
		int u = U[e], v = V[e];
		dg[u] += D[e], dg[v] += D[e];
		if (u) A[u-1][u-1] += D[e];
		if (v) A[v-1][v-1] += D[e];
		if (u && v) A[u-1][v-1] -= D[e], A[v-1][u-1] -= D[e];
	}
	dg[0] = 1;
	auto x = conjGrad([&](const vector<T>& p, vector<T>& y) {
		rep(i,0,n) y[i] = dg[i] * p[i];
		rep(e,0,sz(U)) y[U[e]] -= D[e] * p[V[e]], y[V[e]] -= D[e] * p[U[e]];
		y[0] = 0;
	}, b, dg);
	assert(x[0] == 0);
	vd bd(b.begin() + 1, b.end()), xd(n - 1);
	assert(solveLinear(A, bd, xd) == n - 1);
	rep(i,1,n) assert(abs((double)x[i] - xd[i-1]) <= 1e-7 * max(1.0, abs(xd[i-1])));
}

int main() {
	rep(it,0,3000) {
		int n = rand() % 12 + 1;
		testDense<double>(n);
		testDense<long double>(n);
		int k = rand() % 40 + 2;
		testLaplacian<double>(k, rand() % (2 * k));
		testLaplacian<long double>(k, rand() % (2 * k));
	}
	{ // zero right-hand side
		vd b(3), d(3, 1);
		auto x = conjGrad([&](const vd& v, vd& y) { y = v; }, b, d);
		assert(x == vd(3));
	}
	cout << "Tests passed!" << endl;
}
