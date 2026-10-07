#include "../utilities/template.h"
#include "../utilities/genTree.h"

#include "../../content/graph/HLD.h"

namespace old {
#include "oldHLD.h"
}
struct bruteforce { // values in nodes
    vector<vector<int>> tree;
    vector<int> vals;
    vector<int> pars;
    int unit = -1e9;
    int f(int a, int b) { return max(a, b); }
    void root(int cur, int p = -1) {
        pars[cur] = p;
        for (auto i: tree[cur]) {
            if (i != p) root(i, cur);
        }
    }
    bruteforce(vector<vector<int>> _tree): tree(_tree), vals(sz(tree)), pars(sz(tree)) {
        root(0);
    }
    bool dfsModify(int cur, int target, int val, int p=-1) {
        if (cur == target) {
            vals[cur] += val;
            return true;
        }
        bool alongPath = false;
        for (auto i: tree[cur]) {
            if (i == p) continue;
            alongPath |= dfsModify(i, target, val, cur);
        }
        if (alongPath) vals[cur] += val;
        return alongPath;
    }
    void modifyPath(int a, int b, int val) {
        dfsModify(a, b, val);
    }

    int dfsQuery(int cur, int target, int p = -1) {
        if (cur == target) {
            return vals[cur];
        }
        int res = unit;
        for (auto i: tree[cur]) {
            if (i == p) continue;
            res = f(res, dfsQuery(i, target, cur));
        }
        if (res != unit) {
            return f(res, vals[cur]);
        }
        return res;
    }
    int queryPath(int a, int b) {
        return dfsQuery(a, b);
    }
    int dfsSubtree(int cur, int p) {
        int res = vals[cur];
        for (auto i: tree[cur]) {
            if (i != p)
                res = f(res, dfsSubtree(i, cur));
        }
        return res;
    }
    int querySubtree(int a) {
        return dfsSubtree(a, pars[a]);
    }
};

// Brute force for both node values (E = 0) and edge values (E = 1, the value
// of edge (par[v], v) lives in v), on arbitrary tree shapes rooted at 0.
template <bool E> struct Brute2 {
    int n; vi par, dep, val, res; // res is reused: operator new never frees
    Brute2(const vector<vi>& adj) : n(sz(adj)), par(n, -1), dep(n), val(n) {
        res.reserve(n);
        vi st = {0};
        while (!st.empty()) {
            int v = st.back(); st.pop_back();
            for (int u : adj[v]) if (u != par[v])
                par[u] = v, dep[u] = dep[v] + 1, st.push_back(u);
        }
    }
    vi& path(int a, int b) { // nodes (or child endpoints of edges) on the path
        res.clear();
        while (a != b) {
            if (dep[a] < dep[b]) swap(a, b);
            res.push_back(a); a = par[a];
        }
        if (!E) res.push_back(a);
        return res;
    }
    void modifyPath(int a, int b, int v) { for (int x : path(a, b)) val[x] += v; }
    int queryPath(int a, int b) {
        int r = -1e9; for (int x : path(a, b)) r = max(r, val[x]);
        return r;
    }
    int querySubtree(int v) {
        int r = -1e9;
        rep(x,0,n) {
            int y = x; bool in = 0;
            while (y != -1) { if (y == v) in = 1; y = par[y]; }
            if (in && !(E && x == v)) r = max(r, val[x]);
        }
        return r;
    }
};

// shape: 0 random, 1 path, 2 star, 3 caterpillar, 4 binary; labels shuffled
// except that shape 1/2 with shuf = 0 keeps 0 as the path end / star center.
vector<vi> genShape(int n, int shape, bool shuf) {
    vi lab(n); iota(all(lab), 0);
    if (shuf) random_shuffle(all(lab));
    vector<vi> adj(n);
    rep(i,1,n) {
        int p;
        if (shape == 0) p = rand() % i;
        else if (shape == 1) p = i - 1;
        else if (shape == 2) p = 0;
        else if (shape == 3) p = i % 2 ? i - 1 : max(i - 2, 0);
        else p = (i - 1) / 2;
        adj[lab[i]].push_back(lab[p]);
        adj[lab[p]].push_back(lab[i]);
    }
    for (auto& a : adj) random_shuffle(all(a));
    return adj;
}

template <bool E> void testShapes(int maxn, int iters, int queries) {
    for (int it = 0; it < iters; it++) {
        int n = rand() % maxn + 1;
        auto adj = genShape(n, rand() % 5, rand() % 4 != 0);
        HLD<E> hld(adj);
        Brute2<E> br(adj);
        hld.tree->set(0, n, 0);
        // structural invariants: light edges on any root path <= log2(n)
        rep(v,0,n) {
            int light = 0;
            for (int x = v; x != 0; x = hld.par[x]) {
                assert(hld.par[x] == br.par[x]);
                light += hld.rt[x] == x;
            }
            assert((1 << light) <= n);
            assert(hld.depth[v] == br.dep[v]);
        }
        for (int q = 0; q < queries; q++) {
            int t = rand() % 3, a = rand() % n, b = rand() % n;
            if (t == 0) {
                int val = rand() % 21 - 10;
                hld.modifyPath(a, b, val); br.modifyPath(a, b, val);
            } else if (t == 1) {
                assert(hld.queryPath(a, b) == br.queryPath(a, b));
            } else {
                assert(hld.querySubtree(a) == br.querySubtree(a));
            }
        }
    }
}

void testAgainstOld(int n, int iters, int queries) {
    for (int trees = 0; trees < iters; trees++) {
        auto graph = genRandomTree(n);
        vector<vector<int>> tree1(n);
        vector<vector<pair<int, int>>> tree2(n);
        for (auto i : graph) {
            tree1[i.first].push_back(i.second);
            tree1[i.second].push_back(i.first);
        }
        for (int i = 0; i < sz(tree1); i++) {
            for (auto j : tree1[i]) {
                tree2[i].push_back({j, 0});
            }
        }
        HLD<false> hld(tree1);
        old::HLD hld2(tree2);
        hld.tree->set(0, n, 0);
        for (int itr = 0; itr < queries; itr++) {
            if (rand() % 2) {
                int node = rand() % n;
                int val = rand() % 10;
                hld2.update(node, val);
                hld.modifyPath(node, node, val - hld.queryPath(node, node));
            } else {
                int a = rand() % n;
                int b = rand() % n;
                assert(hld.queryPath(a, b) == hld2.query2(a, b).first);
            }
        }
    }
}
void testAgainstBrute(int n, int iters, int queries) {
    for (int trees = 0; trees < iters; trees++) {
        auto graph = genRandomTree(n);
        vector<vector<int>> tree1(n);
        for (auto i : graph) {
            tree1[i.first].push_back(i.second);
            tree1[i.second].push_back(i.first);
        }
        HLD<false> hld(tree1);
        bruteforce hld2(tree1);
        hld.tree->set(0, n, 0);
        for (int itr = 0; itr < queries; itr++) {
            int rng = rand() % 3;
            if (rng == 0) {
                int a = rand() % n;
                int b = rand() % n;
                int val = rand() % 10;
                hld.modifyPath(a, b, val);
                hld2.modifyPath(a, b, val);
            } else if (rng == 1){
                int a = rand() % n;
                int b = rand() % n;
                hld.queryPath(a, b);
                hld2.queryPath(a, b);
                assert(hld.queryPath(a, b) == hld2.queryPath(a, b));
            } else if (rng == 2) {
                int a = rand() % n;
                assert(hld.querySubtree(a) == hld2.querySubtree(a));
            }
        }
    }

}
int main() {
    srand(2);
#ifdef DEEP_PATH // dfsSz/dfsHld recurse n deep: segfaults with an 8 MB stack
    {
        int n = 200000;
        vector<vi> adj(n);
        rep(i,1,n) adj[i-1].push_back(i), adj[i].push_back(i-1);
        HLD<false> hld(adj);
        hld.tree->set(0, n, 0);
        hld.modifyPath(0, n-1, 3);
        assert(hld.queryPath(5, n-5) == 3);
        cout<<"Tests passed!"<<endl;
        return 0; // the bump allocator cannot fit the other tests as well
    }
#endif
    testAgainstBrute(5, 1000, 10000);
    testAgainstBrute(1000, 100, 100);
    testAgainstOld(5, 1000, 100);
    testAgainstOld(10000, 100, 1000);
    testShapes<false>(1, 50, 20);
    testShapes<true>(1, 50, 20);
    testShapes<false>(2, 500, 30);
    testShapes<true>(2, 500, 30);
    testShapes<false>(6, 20000, 60);
    testShapes<true>(6, 20000, 60);
    testShapes<false>(30, 3000, 200);
    testShapes<true>(30, 3000, 200);
    testShapes<false>(300, 100, 1000);
    testShapes<true>(300, 100, 1000);
    cout<<"Tests passed!"<<endl;
    return 0;
}
