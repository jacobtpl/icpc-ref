#include "../utilities/template.h"

#include "../../content/graph/GlobalMinCut.h"

mt19937 rng(31337);
int ri(int a, int b) { return (int)(rng() % (unsigned)(b - a + 1)) + a; }

void check(const vector<vi>& mat) {
	int n = sz(mat);
	auto pa = globalMinCut(mat);
	if (n <= 1) { assert(pa.first == INT_MAX && pa.second.empty()); return; }
	// returned side is a proper non-empty subset whose cut weight matches
	ll side = 0;
	assert(!pa.second.empty() && sz(pa.second) < n);
	for (int x : pa.second) {
		assert(0 <= x && x < n && !(side >> x & 1));
		side |= 1LL << x;
	}
	ll cutw = 0;
	rep(i,0,n) rep(j,0,n) if ((side >> i & 1) && !(side >> j & 1)) cutw += mat[i][j];
	assert(cutw == pa.first);
	if (n > 16) return;
	// brute force over all bipartitions (vertex 0 fixed on one side)
	ll best = LLONG_MAX;
	vector<ll> row(n);
	for (int m = 1; m < (1 << n); m += 2) if (m != (1 << n) - 1) {
		ll c = 0;
		rep(i,0,n) if (m >> i & 1) rep(j,0,n) if (!(m >> j & 1)) c += mat[i][j];
		best = min(best, c);
	}
	assert(best == pa.first);
}

vector<vi> gen(int n, int maxw, bool diag) {
	vector<vi> mat(n, vi(n));
	int dens = ri(0, 3) ? ri(0, 100) : 100;
	rep(i,0,n) rep(j,0,i) if (ri(0, 99) < dens)
		mat[i][j] = mat[j][i] = ri(0, maxw);
	if (diag) rep(i,0,n) mat[i][i] = ri(0, maxw); // self-loops must not matter
	return mat;
}

int main() {
	check({}); check({{0}}); check({{5}});
	check({{0, 0}, {0, 0}}); check({{0, 3}, {3, 0}});
	rep(it,0,100000) check(gen(ri(2, 5), 3, ri(0, 1)));
	rep(it,0,20000) check(gen(ri(2, 8), ri(1, 10), ri(0, 1)));
	rep(it,0,2000) check(gen(ri(2, 11), 1, 0));
	rep(it,0,300) check(gen(ri(12, 14), 1000, 1));
	rep(it,0,300) check(gen(ri(2, 10), INT_MAX / 100, 0)); // total weight < 2^31
	rep(it,0,200) check(gen(ri(17, 60), 1000, 1)); // cut consistency only
	cout<<"Tests passed!"<<endl;
}
