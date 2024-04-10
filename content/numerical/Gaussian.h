/**
 * Author: Benq
 * Description: Converts A to rref using Gaussian elim.
 * Time: O(N^3)
 */
using T = double; // or mint
const double EPS = 1e-9; // adjust?
int getRow(vector<vector<double>>& m, int R, int i, int nex) {
	pair<double,int> bes{0,-1}; // find row with max abs value
	rep(j,nex,R) ckmax(bes,{abs(m[j][i]),j}); 
	return bes.first < EPS ? -1 : bes.second; }
int getRow(vector<vector<mint>>& m, int R, int i, int nex) {
	rep(j,nex,R) if (m[j][i] != 0) return j;
	return -1; }
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