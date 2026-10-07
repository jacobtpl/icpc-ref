#include "../utilities/template.h"

#include "../../content/graph/BinaryLifting.h"

mt19937 rng(2024);
int ri(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }

void check(const vi& par0, int root, int queries) {
	int n = sz(par0);
	vi P = par0, depth(n, -1);
	depth[root] = 0;
	rep(i,0,n) { // depths by walking up
		vi path;
		int x = i;
		while (depth[x] == -1) path.push_back(x), x = P[x];
		for (int j = sz(path); j--;) depth[path[j]] = depth[P[path[j]]] + 1;
	}
	vector<vi> tbl = treeJump(P);
	assert(P == par0);
	rep(it,0,queries) {
		int a = ri(0, n-1), b = ri(0, n-1);
		// jumps: any step count below 2^sz(tbl), including past the root
		int steps = ri(0, 3) ? ri(0, depth[a] + 2) : ri(0, (1 << sz(tbl)) - 1);
		int x = a;
		rep(i,0,min(steps, depth[a])) x = P[x];
		assert(jmp(tbl, a, steps) == x);
		// lca by walking
		int u = a, v = b;
		while (u != v) {
			if (depth[u] < depth[v]) swap(u, v);
			u = P[u];
		}
		assert(lca(tbl, depth, a, b) == u);
		assert(lca(tbl, depth, a, a) == a);
		assert(lca(tbl, depth, a, root) == root);
	}
}

int main() {
	rep(it,0,50000) {
		int n = ri(1, 12);
		// random tree with random labels, random root
		vi perm(n);
		iota(all(perm), 0);
		shuffle(all(perm), rng);
		vi P(n);
		P[perm[0]] = perm[0];
		int kind = ri(0, 3);
		rep(i,1,n) {
			int p = kind == 0 ? i-1 : kind == 1 ? 0 : ri(max(0, i-3), i-1);
			if (kind == 3) p = ri(0, i-1);
			P[perm[i]] = perm[p];
		}
		check(P, perm[0], 30);
	}
	for (int n : {1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 1000, 1024, 1025}) {
		vi P(n); // path: worst-case depth, exact powers of two
		rep(i,0,n) P[i] = max(i-1, 0);
		check(P, 0, 2000);
		rep(i,0,n) P[i] = min(i+1, n-1);
		check(P, n-1, 2000);
		rep(i,0,n) P[i] = i ? ri(max(0, i-5), i-1) : 0;
		check(P, 0, 2000);
	}
	cout << "Tests passed!" << endl;
}
