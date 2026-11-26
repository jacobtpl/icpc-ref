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
 * Time: O(n m min(n, m)); ker adds O(rank^2 (m - rank))
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
	// Pivot weight; types without abs (mint) take any nonzero pivot.
	auto mag = [](const T& v) {
		if constexpr (requires { abs(v) > abs(v); }) return abs(v);
		else return v != 0;
	};
	auto tol = mag(T()); // weights <= tol count as zero
	if constexpr (is_floating_point<T>::value) tol = T(eps);
	rep(i,0,n) {
		auto v = tol, bv = tol;
		rep(r,i,n) rep(c,i,m)
			if ((v = mag(A[r][c])) > bv)
				br = r, bc = c, bv = v;
		if (!(bv > tol)) {
			rep(j,i,n) if (mag(b[j]) > tol) return -1;
			break;
		}
		swap(A[i], A[br]);
		swap(b[i], b[br]);
		swap(col[i], col[bc]);
		rep(j,0,n) swap(A[j][i], A[j][bc]);
		T inv = 1 / A[i][i];
		rep(j,i+1,n) {
			T fac = A[j][i] * inv;
			b[j] -= fac * b[i];
			rep(k,i+1,m) A[j][k] -= fac * A[i][k];
		}
		rank++;
	}

	x.assign(m, 0);
	if (ker) {
		ker->assign(m - rank, vector<T>(m));
		rep(k,rank,m) (*ker)[k - rank][col[k]] = 1;
	}
	for (int i = rank; i--;) {
		T inv = 1 / A[i][i];
		b[i] *= inv;
		x[col[i]] = b[i];
		rep(j,0,i) b[j] -= A[j][i] * b[i];
		if (ker) { // the same back substitution for each free column
			rep(k,rank,m) if (mag(A[i][k] *= inv) > tol)
				(*ker)[k - rank][col[i]] = -A[i][k];
			rep(j,0,i) rep(k,rank,m) A[j][k] -= A[j][i] * A[i][k];
		}
	}
	return rank; // (multiple solutions if rank < m)
}
