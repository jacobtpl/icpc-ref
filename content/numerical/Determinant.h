/**
 * Author: Simon Lindholm
 * Date: 2016-09-06
 * License: CC0
 * Source: folklore
 * Description: Determinant of a square matrix over double, mint (prime modulus)
 *  or Fraction. Destroys the matrix. For a non-prime modulus use IntDeterminant.h.
 * Time: O(N^3)
 * Status: tested against permutation expansion for mint and Fraction
 */
#pragma once

template<class T>
T det(vector<vector<T>>& a) {
	int n = sz(a); T res = 1;
	rep(i,0,n) {
		int b = i;
		rep(j,i+1,n) {
			if constexpr (is_floating_point<T>::value) {
				if (abs(a[j][i]) > abs(a[b][i])) b = j;
			} else if (a[b][i] == 0) b = j;
		}
		if (i != b) swap(a[i], a[b]), res = -res;
		res *= a[i][i];
		if (res == 0) return 0;
		T inv = 1 / a[i][i];
		rep(j,i+1,n) {
			T v = a[j][i] * inv;
			if (v != 0) rep(k,i+1,n) a[j][k] -= v * a[i][k];
		}
	}
	return res;
}
