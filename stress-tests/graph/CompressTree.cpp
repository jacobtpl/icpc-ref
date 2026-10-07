#include "../utilities/template.h"

#include "../../content/graph/CompressTree.h"

mt19937 rng(777);
int ri(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }

void check(int n, const vector<vi>& adjIn, int iters) {
	vector<vi> adj = adjIn;
	vi par(n, -1), dep(n), order = {0};
	rep(i,0,sz(order)) {
		int v = order[i];
		for (int y : adj[v]) if (y != par[v])
			par[y] = v, dep[y] = dep[v] + 1, order.push_back(y);
	}
	assert(sz(order) == n);
	auto slowLca = [&](int a, int b) {
		while (a != b) {
			if (dep[a] < dep[b]) swap(a, b);
			a = par[a];
		}
		return a;
	};
	LCA lca(adj);
	rep(it,0,iters) {
		int k = ri(1, ri(0, 3) ? min(n, 6) : 2 * n);
		vi subset(k);
		for (int& x : subset) x = ri(0, n-1); // duplicates allowed
		if (ri(0, 1)) {
			sort(all(subset));
			subset.erase(unique(all(subset)), subset.end());
		}
		// expected node set: S closed under pairwise LCA
		set<int> want(all(subset));
		for (int a : subset) for (int b : subset) want.insert(slowLca(a, b));
		vpi ret = compressTree(lca, subset);
		assert(sz(ret) == sz(want));
		set<int> got;
		for (auto& pr : ret) got.insert(pr.second);
		assert(got == want);
		if (sz(set<int>(all(subset))) > 1) assert(sz(ret) <= 2 * sz(subset) - 1);
		assert(ret[0].first == 0); // root points to itself
		rep(i,1,sz(ret)) {
			// parent comes earlier and is the nearest proper ancestor in the set
			assert(0 <= ret[i].first && ret[i].first < i);
			int x = par[ret[i].second];
			while (!want.count(x)) x = par[x];
			assert(ret[ret[i].first].second == x);
		}
	}
}

int main() {
	rep(it,0,30000) {
		int n = ri(1, 12), kind = ri(0, 3);
		vi perm(n);
		iota(all(perm), 0);
		shuffle(perm.begin() + 1, perm.end(), rng); // LCA.h roots at 0
		vector<vi> adj(n);
		rep(i,1,n) {
			int p = kind == 0 ? i-1 : kind == 1 ? 0 : ri(max(0, i-3), i-1);
			if (kind == 3) p = ri(0, i-1);
			adj[perm[i]].push_back(perm[p]);
			adj[perm[p]].push_back(perm[i]);
		}
		for (auto& v : adj) shuffle(all(v), rng);
		check(n, adj, 20);
	}
	rep(it,0,30) {
		int n = ri(100, 300);
		vector<vi> adj(n);
		rep(i,1,n) {
			int p = ri(max(0, i - (it % 3 == 0 ? 1 : it % 3 == 1 ? 5 : n)), i-1);
			adj[i].push_back(p);
			adj[p].push_back(i);
		}
		check(n, adj, 40);
	}
	cout << "Tests passed!" << endl;
}
