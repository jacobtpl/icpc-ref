/**
 * Author: Ulf Lundstrom, Simon Lindholm
 * Date: 2009-08-15
 * License: CC0
 * Source: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
 * Description: x = tridiagonal(d, p, q, b) solves the n x n system
 *   d[0] x[0] + p[0] x[1] = b[0]
 *   q[i-1] x[i-1] + d[i] x[i] + p[i] x[i+1] = b[i]   (0 < i < n-1)
 *   q[n-2] x[n-2] + d[n-1] x[n-1] = b[n-1]
 * i.e. d is the diagonal, p the superdiagonal and q the subdiagonal.
 * To solve a_i = b_i a_(i-1) + c_i a_(i+1) + d_i for 1 <= i <= n with
 * a_0 and a_(n+1) known, use
 *   a = tridiagonal({1,-1,...,-1,1}, {0,c_1,...,c_n},
 *                   {b_1,...,b_n,0}, {a_0,-d_1,...,-d_n,a_(n+1)}).
 * Fails if the solution is not unique.
 * If |d_i| > |p_i| + |q_(i-1)| for all i, or |d_i| > |p_(i-1)| + |q_i|, or the matrix is positive definite,
 * the algorithm is numerically stable and neither tr nor the check for diag[i] == 0 is needed.
 * Time: O(N)
 * Status: Brute-force tested mod 5 and 7 and stress-tested for real matrices obeying the criteria above.
 */
#pragma once

typedef double T;
template<class T>
vector<T> tridiagonal(vector<T> diag, const vector<T>& super,
		const vector<T>& sub, vector<T> b) {
	int n = sz(b); vi tr(n);
	rep(i,0,n-1) {
		bool zero = diag[i] == 0;
		if constexpr (is_floating_point<T>::value)
			zero = abs(diag[i]) < 1e-9 * abs(super[i]);
		if (zero) {
			b[i+1] -= b[i] * diag[i+1] / super[i];
			if (i+2 < n) b[i+2] -= b[i] * sub[i+1] / super[i];
			diag[i+1] = sub[i]; tr[++i] = 1;
		} else {
			diag[i+1] -= super[i]*sub[i]/diag[i];
			b[i+1] -= b[i]*sub[i]/diag[i];
		}
	}
	for (int i = n; i--;) {
		if (tr[i]) {
			swap(b[i], b[i-1]);
			diag[i-1] = diag[i];
			b[i] /= super[i-1];
		} else {
			b[i] /= diag[i];
			if (i) b[i-1] -= b[i]*super[i-1];
		}
	}
	return b;
}
