#include "../utilities/template.h"
#include "../utilities/genTree.h"

#include "../../content/graph/LCA.h"
#include "../../content/graph/BinaryLifting.h"
#include "../../content/data-structures/RMQ.h"

namespace old {
typedef vector<pii> vpi;
typedef vector<vpi> graph;

struct LCA {
    vi time;
    vector<ll> dist;
    RMQ<pii> rmq;

    LCA(graph& C) : time(sz(C), -99), dist(sz(C)), rmq(dfs(C)) {}

    vpi dfs(graph& C) {
        vector<tuple<int, int, int, ll>> q(1);
        vpi ret;
        int T = 0, v, p, d; ll di;
        while (!q.empty()) {
            tie(v, p, d, di) = q.back();
            q.pop_back();
            if (d) ret.emplace_back(d, p);
            time[v] = T++;
            dist[v] = di;
            for(auto &e: C[v]) if (e.first != p)
                q.emplace_back(e.first, v, d+1, di + e.second);
        }
        return ret;
    }

    int query(int a, int b) {
        if (a == b) return a;
        a = time[a], b = time[b];
        return rmq.query(min(a, b), max(a, b)).second;
    }
    ll distance(int a, int b) {
        int lca = query(a, b);
        return dist[a] + dist[b] - 2 * dist[lca];
    }
};
}


void getPars(vector<vi> &tree, int cur, int p, int d, vector<int> &par, vector<int> &depth) {
    par[cur] = p;
    depth[cur] = d;
    for(auto i: tree[cur]) if (i != p) {
        getPars(tree, i, cur, d+1, par, depth);
    }
}
void test_n(int n, int num) {
    for (int out=0; out<num; out++) {
        auto graph = genRandomTree(n);
        vector<vi> tree(n);
        vector<vector<pair<int, int>>> oldTree(n);
        for (auto i: graph) {
            tree[i.first].push_back(i.second);
            tree[i.second].push_back(i.first);
            oldTree[i.first].push_back({i.second, 1});
            oldTree[i.second].push_back({i.first, 1});
        }
        vector<int> par(n), depth(n);
        getPars(tree, 0, 0, 0, par, depth);
        vector<vi> tbl = treeJump(par);
        LCA new_lca(tree);
        old::LCA old_lca(oldTree);
        for (int i=0; i<100; i++) {
            int a = rand()%n, b = rand()%n;
            int binLca = lca(tbl, depth, a, b);
            int newLca = new_lca.lca(a,b);
            int oldLca = old_lca.query(a,b);
            assert(oldLca == newLca);
            assert(binLca == newLca);
        }
    }
}

// Oracle: walk up parent pointers. `directed` passes child-only lists.
void test_naive(int n, bool directed, int shape) {
	vi par(n, -1), depth(n), perm(n);
	iota(all(perm), 0);
	random_shuffle(perm.begin() + 1, perm.end()); // 0 stays root
	vector<vi> tree(n);
	rep(i,1,n) {
		int p = shape == 0 ? rand() % i : shape == 1 ? i - 1 :
			shape == 2 ? 0 : (i - 1) / 2;
		int a = perm[i], b = perm[p];
		par[a] = b;
		tree[b].push_back(a);
		if (!directed) tree[a].push_back(b);
	}
	if (rand() % 2) for (auto& v : tree) random_shuffle(all(v));
	rep(i,1,n) depth[perm[i]] = depth[par[perm[i]]] + 1;
	LCA l(tree);
	auto naive = [&](int a, int b) {
		while (a != b) {
			if (depth[a] < depth[b]) swap(a, b);
			a = par[a];
		}
		return a;
	};
	if (n <= 40) { rep(a,0,n) rep(b,0,n) assert(l.lca(a, b) == naive(a, b)); }
	else rep(it,0,300) {
		int a = rand() % n, b = rand() % n;
		assert(l.lca(a, b) == naive(a, b));
	}
}

signed main() {
    srand(5);
    rep(n,1,13) rep(it,0,1500) test_naive(n, it & 1, 0);
    rep(it,0,2000) test_naive(rand() % 40 + 1, it & 1, it % 4);
    rep(it,0,40) test_naive(rand() % 3000 + 1, it & 1, it % 4);
    test_n(10, 1000);
    test_n(100, 100);
    test_n(1000, 10);
#ifdef LCA_DEEP
    // Opt-in: the recursive dfs needs ~32 bytes of stack per level, so a
    // path with >= ~3e5 nodes overflows a default 8 MB stack (run with
    // `ulimit -s 8192`; run-all.sh raises the limit to 512 MB).
    {
        int n = 1000000;
        vector<vi> tree(n);
        rep(i,1,n) tree[i-1].push_back(i), tree[i].push_back(i-1);
        LCA l(tree);
        rep(i,0,1000) {
            int a = rand() % n, b = rand() % n;
            assert(l.lca(a, b) == min(a, b));
        }
    }
#endif
    cout<<"Tests passed!"<<endl;
}

