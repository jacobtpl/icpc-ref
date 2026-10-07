#include "../utilities/template.h"
#include "../utilities/utils.h"

#include "../../content/graph/GlobalMinCut.h"
#include "../../content/graph/GomoryHu.h"
// Independent max-flow oracle (adjacency-matrix Dinic). content/graph/Dinic.h
// cannot be used here: it defines a global `struct Edge` that clashes with
// GomoryHu.h's `Edge` typedef and has a different interface.
struct Dinic {
    int n; vector<vector<ll>> c; vi lvl, ptr;
    Dinic(int n) : n(n), c(n, vector<ll>(n)), lvl(n), ptr(n) {}
    void addEdge(int a, int b, ll cap, ll rcap = 0) {
        if (a != b) c[a][b] += cap, c[b][a] += rcap;
    }
    ll dfs(int v, int t, ll f) {
        if (v == t || !f) return f;
        for (int& u = ptr[v]; u < n; u++) if (lvl[u] == lvl[v] + 1 && c[v][u])
            if (ll p = dfs(u, t, min(f, c[v][u]))) {
                c[v][u] -= p, c[u][v] += p;
                return p;
            }
        return 0;
    }
    ll calc(int s, int t) {
        ll flow = 0;
        for (;;) {
            fill(all(lvl), -1); fill(all(ptr), 0);
            vi q = {s}; lvl[s] = 0;
            rep(i,0,sz(q)) rep(u,0,n) if (lvl[u] < 0 && c[q[i]][u])
                lvl[u] = lvl[q[i]] + 1, q.push_back(u);
            if (lvl[t] < 0) return flow;
            while (ll p = dfs(s, t, LLONG_MAX)) flow += p;
        }
    }
};


void test(int N, int mxFlow, int iters) {
    for (int it = 0; it < iters; it++) {
        int n = rand()%N+1;
        int m = rand()%(N*N);
        vector<array<ll, 3>> edges;
        vector<vi> mat(n, vi(n));
        rep(it,0,m) {
            int i = rand() % n;
            int j = rand() % n;
            if (i == j) continue;
            int w = rand() % mxFlow;
            edges.push_back({i, j, w});
            mat[i][j] += w;
            mat[j][i] += w;
        }
        auto calc = [&](int s, int t) {
            Dinic flow(n);
            for (auto e : edges) {
                flow.addEdge((int)e[0], (int)e[1], e[2], e[2]);
            }
            return flow.calc(s, t);
        };
        vector<Edge> gomoryHuTree = gomoryHu(n, edges);
        vector<vector<array<int, 2>>> adj(n);
        for (auto e : gomoryHuTree) {
            adj[e[0]].push_back({(int)e[1], (int)e[2]});
            adj[e[1]].push_back({(int)e[0], (int)e[2]});
        }
        auto dfs = make_y_combinator([&](auto dfs, int start, int cur, int p, int mn) -> void {
            if (start != cur) {
                assert(mn == calc(start, cur));
            }
            for (auto i : adj[cur]) {
                if (i[0] != p)
                    dfs(start, i[0], cur, min(mn, i[1]));
            }
        });
        dfs(0, 0, -1, INT_MAX);

        // Check that the lightest edge agrees with GlobalMinCut.
        if (n >= 2) {
            ll minCut = LLONG_MAX;
            for (auto e : gomoryHuTree) {
                minCut = min(minCut, e[2]);
            }
            auto mat2 = mat;
            auto pa = globalMinCut(mat2);
            assert(pa.first == minCut);
            vi inCut(n);
            assert(sz(pa.second) != 0);
            assert(sz(pa.second) != n);
            for (int x : pa.second) {
                assert(0 <= x && x < n);
                assert(!inCut[x]);
                inCut[x] = 1;
            }
            int cutw = 0;
            rep(i,0,n) rep(j,0,n) if (inCut[i] && !inCut[j]) {
                cutw += mat[i][j];
            }
            assert(pa.first == cutw);
        }
    }
}
// All-pairs check against Dinic, including self-loops, multi-edges, zero
// capacities, disconnected graphs and capacities needing 64 bits.
void testAllPairs(int N, ll mxFlow, int iters) {
    for (int it = 0; it < iters; it++) {
        int n = rand()%N+1;
        int m = rand()%(n*n+1);
        if (rand()%3 == 0) m = rand()%(n+1); // sparse / disconnected
        vector<array<ll, 3>> edges;
        rep(it,0,m) {
            int i = rand() % n, j = rand() % n; // i == j allowed
            ll w = (ll)((((unsigned long long)rand() << 31) ^ rand()) % mxFlow);
            edges.push_back({i, j, w});
        }
        vector<Edge> tree = gomoryHu(n, edges);
        assert(sz(tree) == n - 1);
        vector<vector<ll>> mn(n, vector<ll>(n, LLONG_MAX));
        vector<vector<pair<int, ll>>> adj(n);
        for (auto e : tree) {
            assert(0 <= e[0] && e[0] < n && 0 <= e[1] && e[1] < n);
            adj[e[0]].push_back({(int)e[1], e[2]});
            adj[e[1]].push_back({(int)e[0], e[2]});
        }
        rep(s,0,n) { // min edge on tree path from s (also checks it is a tree)
            vi seen(n), st = {s}; seen[s] = 1; int cnt = 0;
            while (!st.empty()) {
                int v = st.back(); st.pop_back(); cnt++;
                for (auto [u, w] : adj[v]) if (!seen[u])
                    seen[u] = 1, mn[s][u] = min(mn[s][v], w), st.push_back(u);
            }
            assert(cnt == n);
        }
        rep(s,0,n) rep(t,0,s) {
            Dinic flow(n);
            for (auto e : edges) flow.addEdge((int)e[0], (int)e[1], e[2], e[2]);
            assert(flow.calc(s, t) == mn[s][t]);
        }
    }
}
signed main() {
    test(25, 5, 200);
    test(100, 1000, 5);
    test(100, 1, 20);
    test(5, 5, 20000);
    testAllPairs(1, 5, 10);
    testAllPairs(4, 3, 20000);
    testAllPairs(7, 4, 5000);
    testAllPairs(12, 100, 500);
    testAllPairs(8, (ll)1e15, 2000);
    testAllPairs(40, 1000000, 10);
    cout<<"Tests passed!"<<endl;
}
