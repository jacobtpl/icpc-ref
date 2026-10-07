#include "../utilities/template.h"
#include "../utilities/random.h"

#include "../../content/graph/TopoSort.h"

int main() {
	rep(it,0,50000) {
		int n = rand() % 20;
		int m = n ? rand() % 30 : 0;
		bool acyclic = randBool();
		vi order(n);
		iota(all(order), 0);
		shuffle_vec(order);
		vector<vi> ed(n);
		rep(i,0,m) {
			int a = rand() % n;
			int b = rand() % n;
			if (acyclic && a >= b) continue;
			ed[order[a]].push_back(order[b]);
		}
		vi ret = topoSort(ed);
		if (acyclic) assert(sz(ret) == n);
		else assert(sz(ret) <= n);
		vi seen(n);
		for (int i : ret) {
			assert(!seen[i]++);
			for (int j : ed[i])
				assert(!seen[j]);
		}
	}
	// exact contract: the result is precisely the set of nodes not reachable from a cycle
	// (self-loops and multi-edges included)
	assert(topoSort({}).empty());
	assert(topoSort({{}}) == vi{0});
	assert(topoSort({{0}}).empty());
	assert(topoSort({{1, 1}, {}}) == (vi{0, 1}));
	rep(it,0,50000) {
		int n = rand() % 9 + 1, m = rand() % 14;
		vector<vi> ed(n);
		vector<vi> reach(n, vi(n));
		rep(i,0,m) {
			int a = rand() % n, b = rand() % n;
			ed[a].push_back(b); reach[a][b] = 1;
		}
		rep(k,0,n) rep(i,0,n) rep(j,0,n) if (reach[i][k] && reach[k][j]) reach[i][j] = 1;
		vi bad(n);
		rep(i,0,n) if (reach[i][i]) { bad[i] = 1; rep(j,0,n) if (reach[i][j]) bad[j] = 1; }
		vi ret = topoSort(ed), pos(n, -1);
		rep(i,0,sz(ret)) { assert(pos[ret[i]] == -1); pos[ret[i]] = i; }
		rep(i,0,n) assert((pos[i] == -1) == bad[i]);
		rep(i,0,n) if (!bad[i]) for (int j : ed[i]) if (!bad[j]) assert(pos[i] < pos[j]);
	}
	cout << "Tests passed!" << endl;
	return 0;
}
