#include "../utilities/template.h"

#include "../../content/data-structures/SubMatrix.h"

mt19937_64 rng(777);
ll rnd(ll a, ll b) { return a + (ll)(rng() % (unsigned long long)(b - a + 1)); }

template<class T>
void run(int R, int C, ll lo, ll hi) {
	vector<vector<T>> v(R, vector<T>(C));
	for (auto& row : v) for (auto& x : row) x = (T)rnd(lo, hi);
	SubMatrix<T> m(v);
	rep(u,0,R+1) rep(d,u,R+1) rep(l,0,C+1) rep(r,l,C+1) {
		T s = 0;
		rep(i,u,d) rep(j,l,r) s += v[i][j];
		assert(m.sum(u, l, d, r) == s);
	}
}

int main() {
	{ // usage example from the header
		vector<vi> matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
		SubMatrix<int> m(matrix);
		assert(m.sum(0, 0, 2, 2) == 12);
	}
	{ // empty matrix and matrices with empty rows
		vector<vi> e;
		SubMatrix<int> m(e);
		assert(m.sum(0, 0, 0, 0) == 0);
		vector<vi> e2(3);
		SubMatrix<int> m2(e2);
		assert(m2.sum(0, 0, 3, 0) == 0);
	}
	rep(it,0,4000) {
		int R = (int)rnd(1, 7), C = (int)rnd(1, 7);
		run<int>(R, C, -1000, 1000);
		run<ll>(R, C, -(ll)1e16, (ll)1e16); // 49 * 1e16 < 9.2e18
		run<int>(R, C, -40000000, 40000000); // 49 * 4e7 < 2^31
		run<ll>(R, C, 0, 1);
	}
	rep(it,0,3) run<ll>(1, 40, -(ll)1e17, (ll)1e17), run<ll>(40, 1, -(ll)1e17, (ll)1e17);
	rep(it,0,2) run<ll>(25, 25, -(ll)1e15, (ll)1e15);
	cout<<"Tests passed!"<<endl;
}
