/**
 * Author: Devin
 * Date: 2026-11-25
 * License: CC0
 * Source: UCup Season 5 Stage 8 H; Hestenes-Stiefel with Jacobi preconditioner
 * Description: Solves A x = b for a symmetric positive definite A given only
 * as a function: mul(x, y) must set y = A x (y is already sized). Use it when
 * A is large and sparse, so mul is O(nonzeros) and Gaussian elimination is
 * too slow. diag = diagonal of A (speeds up convergence; pass all ones if
 * unknown). Stops when |b - A x| <= eps * |b| or after maxIt iterations.
 * Graph Laplacian (electrical flow, min sum of flow^2 / w): fix one node to 0
 * so A is definite: its b and diag entries are ignored if mul keeps y[0] = 0
 * and b[0] = 0, diag[0] = 1.
 * Use long double on x86 if the system is ill-conditioned (long paths).
 * Time: O(iterations * mul). Iterations <= n in exact arithmetic, about
 * sqrt(condition number) in practice; can exceed n due to rounding.
 * Status: stress-tested
 */
#pragma once

template<class T, class F>
vector<T> conjGrad(F mul, const vector<T>& b, const vector<T>& diag,
		T eps = 1e-12, int maxIt = -1) {
	int n = sz(b);
	if (maxIt < 0) maxIt = 20 * n + 100;
	vector<T> x(n), r = b, z(n), p(n), q(n);
	T rz = 0, bb = 0;
	rep(i,0,n) p[i] = r[i] / diag[i], rz += r[i] * p[i], bb += b[i] * b[i];
	while (maxIt-- && bb > 0) {
		mul(p, q);
		T pq = 0, rr = 0, nz = 0;
		rep(i,0,n) pq += p[i] * q[i];
		T a = rz / pq;
		rep(i,0,n) x[i] += a * p[i], r[i] -= a * q[i], rr += r[i] * r[i];
		if (rr <= eps * eps * bb) break;
		rep(i,0,n) z[i] = r[i] / diag[i], nz += r[i] * z[i];
		T be = nz / rz; rz = nz;
		rep(i,0,n) p[i] = z[i] + be * p[i];
	}
	return x;
}
