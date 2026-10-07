#include "../utilities/template.h"

#include "../../content/graph/Centroid.h"

mt19937 rng(4242);
int ri(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }

vector<vi> genTree(int n, int kind) {
	vi perm(n);
	iota(all(perm), 0);
	shuffle(all(perm), rng);
	vector<vi> adj(n);
	rep(i,1,n) {
		int p = kind == 0 ? i-1 : kind == 1 ? 0 : kind == 2 ? ri(max(0, i-3), i-1)
			: kind == 3 ? ri(0, i-1) : kind == 4 ? (i-1) / 2 : (i % 2 ? i-1 : max(0, i-2));
		adj[perm[i]].push_back(perm[p]);
		adj[perm[p]].push_back(perm[i]);
	}
	for (auto& v : adj) shuffle(all(v), rng);
	return adj;
}

// Returns the depth of the centroid tree after verifying it from scratch.
int check(const vector<vi>& adj) {
	int n = sz(adj);
	Centroid C(adj);
	assert(C.N == n && sz(C.par) == n && sz(C.links) == n);
	assert(C.par[C.root] == pii(-1, -1));
	rep(v,0,n) {
		assert(C.rem[v]);
		rep(i,0,sz(C.links[v])) assert(C.par[C.links[v][i]] == pii(v, i));
		if (v != C.root) {
			auto [p, i] = C.par[v];
			assert(0 <= p && p < n && 0 <= i && i < sz(C.links[p]) && C.links[p][i] == v);
		}
	}
	// replay the decomposition top-down
	vi gone(n), stamp(n, -1);
	int maxDepth = 0, seen = 0, T = 0;
	vector<pii> todo = {{C.root, 1}};
	// nodes of the component of `from` among non-removed nodes
	auto collect = [&](int from) {
		vi q = {from};
		stamp[from] = ++T;
		rep(i,0,sz(q)) for (int y : adj[q[i]]) if (!gone[y] && stamp[y] != T)
			stamp[y] = T, q.push_back(y);
		return q;
	};
	while (!todo.empty()) {
		auto [c, d] = todo.back();
		todo.pop_back();
		maxDepth = max(maxDepth, d);
		seen++;
		assert(!gone[c]);
		int total = sz(collect(c));
		if (c == C.root) assert(total == n);
		gone[c] = 1;
		// each remaining neighbour's component: at most half, exactly one child centroid in it
		vi kids;
		for (int x : adj[c]) if (!gone[x]) {
			vi q = collect(x);
			assert(2 * sz(q) <= total);
			int cnt = 0;
			for (int y : q) for (int k : C.links[c]) if (k == y) cnt++, kids.push_back(k);
			assert(cnt == 1);
		}
		assert(sz(kids) == sz(C.links[c]));
		for (int k : C.links[c]) todo.push_back({k, d + 1});
	}
	assert(seen == n);
	return maxDepth;
}

int main() {
	rep(it,0,60000) {
		int n = ri(1, 14);
		int d = check(genTree(n, ri(0, 5)));
		assert((1 << (d - 1)) <= n);
	}
	rep(it,0,60) {
		int n = ri(200, 1500);
		int d = check(genTree(n, it % 6));
		assert((1 << (d - 1)) <= n);
	}
	cout << "Tests passed!" << endl;
}
