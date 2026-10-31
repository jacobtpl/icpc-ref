/**
 * Author: Benq
 * Description: Converts A to rref using Gaussian elim.
 * Time: O(N^3)
 */
using T = double; // or mint
const double EPS = 1e-9; // adjust?
template<class T>
int getRow(vector<vector<T>>& m, int R, int i, int nex) {
	if constexpr (is_floating_point<T>::value) {
		pair<T,int> best{0,-1};
		rep(j,nex,R) best = max(best, {abs(m[j][i]), j});
		return best.first < EPS ? -1 : best.second;
	} else {
		rep(j,nex,R) if (m[j][i] != 0) return j;
		return -1;
	}
}
template<class T>
pair<T,int> gauss(vector<vector<T>>& m) { // convert to rref
	if (!sz(m)) return {1,0};
	int R = sz(m), C = sz(m[0]), rank = 0, nex = 0;
	T prod = 1; // determinant
	rep(i,0,C) {
		int row = getRow(m,R,i,nex);
		if (row == -1) { prod = 0; continue; }
		if (row != nex) prod *= -1, swap(m[row],m[nex]);
		prod *= m[nex][i]; rank++;
		T x = 1/m[nex][i]; rep(k,i,C) m[nex][k] *= x;
		rep(j,0,R) if (j != nex) {
			T v = m[j][i]; if (v == 0) continue;
			rep(k,i,C) m[j][k] -= v*m[nex][k];
		}
		nex++;
	}
	return {prod,rank};
}