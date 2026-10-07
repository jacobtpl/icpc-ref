/**
 * Author: stanford
 * Date: idk
 * License: CC0
 * Source: stanford-notebook
 * Description: Dinic's without scaling. Self-loops are ignored; requires
 * $S \ne T$ (otherwise MaxFlow never returns). The $i$-th added non-loop edge
 * is E[2*i], with its flow in E[2*i].flow. After MaxFlow, $d[v] \le N$ iff
 * $v$ is on the source side of a min cut. Recursion depth is up to $V$.
 * Time: $O(V^2 E)$, $O(\sqrt{V}E)$ for bipartite matching
 * Usage: Dinic d(n); d.AddEdge(u, v, cap); ll f = d.MaxFlow(s, t);
 * Status: stress-tested
 */

struct Dinic {
  struct Edge {
    int u, v;
    ll cap, flow;
    Edge() {}
    Edge(int u, int v, ll cap): u(u), v(v), cap(cap), flow(0) {}
  };
  int N;
  vector<Edge> E;
  vector<vector<int>> g;
  vector<int> d, pt;
  Dinic(int N): N(N), E(0), g(N), d(N), pt(N) {}
  void AddEdge(int u, int v, ll cap) {
    if (u != v) {
      E.emplace_back(u, v, cap);
      g[u].emplace_back(E.size() - 1);
      E.emplace_back(v, u, 0);
      g[v].emplace_back(E.size() - 1);
    }
  }
  bool BFS(int S, int T) {
    queue<int> q({S});
    fill(d.begin(), d.end(), N + 1);
    d[S] = 0;
    while(!q.empty()) {
      int u = q.front(); q.pop();
      if (u == T) break;
      for (int k: g[u]) {
        Edge &e = E[k];
        if (e.flow < e.cap && d[e.v] > d[e.u] + 1) {
          d[e.v] = d[e.u] + 1;
          q.emplace(e.v);
        }
      }
    }
    return d[T] != N + 1;
  }
  ll DFS(int u, int T, ll flow = -1) {
    if (u == T || flow == 0) return flow;
    for (int &i = pt[u]; i < g[u].size(); ++i) {
      Edge &e = E[g[u][i]];
      Edge &oe = E[g[u][i]^1];
      if (d[e.v] == d[e.u] + 1) {
        ll amt = e.cap - e.flow;
        if (flow != -1 && amt > flow) amt = flow;
        if (ll pushed = DFS(e.v, T, amt)) {
          e.flow += pushed;
          oe.flow -= pushed;
          return pushed;
        }
      }
    }
    return 0;
  }
  ll MaxFlow(int S, int T) {
    ll total = 0;
    while (BFS(S, T)) {
      fill(pt.begin(), pt.end(), 0);
      while (ll flow = DFS(S, T))
        total += flow;
    }
    return total;
  }
};