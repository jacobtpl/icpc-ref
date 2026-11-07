/**
 * Author: Simon Lindholm
 * Date: 2016-09-06
 * License: CC0
 * Source: folklore
 * Description: Calculates determinant of a matrix. Destroys the matrix.
 * Time: $O(N^3)$
 * Status: somewhat tested
 */
#pragma once

template<class T>
T det(vector<vector<T>>& a) {
	int n = sz(a); T res = 1;
	rep(i,0,n) {
		int b = i;
		rep(j,i+1,n) if (abs(a[j][i]) > abs(a[b][i])) b = j;
		if (i != b) swap(a[i], a[b]), res *= -1;
		res *= a[i][i];
		if (res == 0) return 0;
		rep(j,i+1,n) {
			T v = a[j][i] / a[i][i];
			if (v != 0) rep(k,i+1,n) a[j][k] -= v * a[i][k];
		}
	}
	return res;
}
