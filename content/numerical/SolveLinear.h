/**
 * Author: Per Austrin, Simon Lindholm
 * Date: 2004-02-08
 * License: CC0
 * Description: Solves A x = b for an n by m matrix A over double, mint or Fraction.
 *  Returns the rank, or -1 if there is no solution. x (which must have size m)
 *  is set to one solution, with the free variables 0.
 *  If ker is given it is set to a basis of the null space (m - rank vectors):
 *  the solutions are exactly x plus any combination of ker, and x[j] is
 *  uniquely determined iff ker[k][j] is zero for every k.
 *  Destroys A and b.
 * Usage: vector<vector<mint>> ker;
 *  int rank = solveLinear(A, b, x, &ker);
 * Time: O(n m min(n, m))
 * Status: brute-force tested mod 2, 3 and 5 for n m <= 4, stress-tested for double and Fraction
 */
#pragma once

const double eps = 1e-12;

template<class T>
int solveLinear(vector<vector<T>>& A, vector<T>& b, vector<T>& x,
		vector<vector<T>>* ker = 0) {
	int n = sz(A), m = sz(x), rank = 0, br = 0, bc = 0;
	if (n) assert(sz(A[0]) == m);
	vi col(m); iota(all(col), 0);
	// Pivot weight: 0 for zero; exact types take the first nonzero.
	auto mag = [](const T& v) -> double {
		if constexpr (is_floating_point<T>::value)
			return abs(v) > eps ? (double)abs(v) : 0;
		else return v != 0;
	};
	rep(i,0,n) {
		double v, bv = 0;
		rep(r,i,n) rep(c,i,m)
			if ((v = mag(A[r][c])) > bv)
				br = r, bc = c, bv = v;
		if (bv == 0) {
			rep(j,i,n) if (mag(b[j]) != 0) return -1;
			break;
		}
		swap(A[i], A[br]);
		swap(b[i], b[br]);
		swap(col[i], col[bc]);
		rep(j,0,n) swap(A[j][i], A[j][bc]);
		T inv = 1 / A[i][i];
		b[i] *= inv;
		rep(k,i,m) A[i][k] *= inv;
		rep(j,0,n) if (j != i) {
			T fac = A[j][i];
			b[j] -= fac * b[i];
			rep(k,i+1,m) A[j][k] -= fac * A[i][k];
		}
		rank++;
	}

	x.assign(m, 0);
	rep(i,0,rank) x[col[i]] = b[i];
	if (ker) {
		ker->assign(m - rank, vector<T>(m));
		rep(j,rank,m) {
			auto& w = (*ker)[j - rank];
			w[col[j]] = 1;
			rep(i,0,rank) w[col[i]] = -A[i][j];
		}
	}
	return rank; // (multiple solutions if rank < m)
}
