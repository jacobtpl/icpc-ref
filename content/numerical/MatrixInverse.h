/**
 * Author: Max Bennedich
 * Date: 2004-02-08
 * Description: Inverts a square matrix A over double, mint (prime modulus)
 *  or Fraction. Returns the rank; the inverse is stored in A unless A is
 *  singular (rank < n), in which case A is destroyed.
 *  For a prime power modulus p^k, start from B = inverse of A mod p and
 *  repeat B = B (2I - A B) mod p^k, doubling k each step.
 * Time: O(n^3)
 * Status: tested by multiplying back for double, mint and Fraction
 */
#pragma once

using T = double; // or mint
template<class T>
int matInv(vector<vector<T>>& A) {
	int n = sz(A); vi col(n);
	vector<vector<T>> tmp(n, vector<T>(n));
	rep(i,0,n) tmp[i][i] = 1, col[i] = i;

	rep(i,0,n) {
		int r = i, c = i;
		rep(j,i,n) rep(k,i,n) {
			// Types without abs (mint) take any nonzero pivot.
			if constexpr (requires { abs(A[j][k]) > abs(A[r][c]); }) {
				if (abs(A[j][k]) > abs(A[r][c])) r = j, c = k;
			} else if (A[r][c] == 0) r = j, c = k;
		}
		bool zero = A[r][c] == 0;
		if constexpr (is_floating_point<T>::value)
			zero = abs(A[r][c]) < 1e-12;
		if (zero) return i;
		A[i].swap(A[r]); tmp[i].swap(tmp[r]);
		rep(j,0,n)
			swap(A[j][i], A[j][c]), swap(tmp[j][i], tmp[j][c]);
		swap(col[i], col[c]);
		T v = 1 / A[i][i];
		rep(j,i+1,n) {
			T f = A[j][i] * v;
			A[j][i] = 0;
			rep(k,i+1,n) A[j][k] -= f*A[i][k];
			rep(k,0,n) tmp[j][k] -= f*tmp[i][k];
		}
		rep(j,i+1,n) A[i][j] *= v;
		rep(j,0,n) tmp[i][j] *= v;
		A[i][i] = 1;
	}
	/// forget A at this point, just eliminate tmp backward
	for (int i = n-1; i > 0; --i) rep(j,0,i) {
		T v = A[j][i];
		rep(k,0,n) tmp[j][k] -= v*tmp[i][k];
	}
	rep(i,0,n) rep(j,0,n) A[col[i]][col[j]] = tmp[i][j];
	return n;
}
