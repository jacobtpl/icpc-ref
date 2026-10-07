#include "../utilities/template.h"

struct Bumpalloc {
	char buf[450 << 20];
	size_t bufp;
	void* alloc(size_t s) {
		assert(s < bufp);
		return (void*)&buf[bufp -= s];
	}
	Bumpalloc() { reset(); }

	template<class T> T* operator=(T&& x) {
		T* r = (T*)alloc(sizeof(T));
		new(r) T(move(x));
		return r;
	}
	void reset() { bufp = sizeof buf; }
} bumpalloc;

// When not testing perf, we don't want to leak memory
#ifndef TEST_PERF
#define new bumpalloc =
#endif
#include "../../content/graph/DirectedMST.h"
#ifndef TEST_PERF
#undef new
#endif

namespace mit {

#define N 110000
#define M 110000
#define inf 2000000000

struct edg {
    int u, v;
    int cost;
} E[M], E_copy[M];

int In[N], ID[N], vis[N], pre[N];

// edges pointed from root.
int Directed_MST(int root, int NV, int NE) {
    for (int i = 0; i < NE; i++)
        E_copy[i] = E[i];
    int ret = 0;
    int u, v;
    while (true) {
        rep(i,0,NV)   In[i] = inf;
        rep(i,0,NE) {
            u = E_copy[i].u;
            v = E_copy[i].v;
            if(E_copy[i].cost < In[v] && u != v) {
                In[v] = E_copy[i].cost;
                pre[v] = u;
            }
        }
        rep(i,0,NV) {
            if(i == root)   continue;
            if(In[i] == inf)    return -1; // no solution
        }

        int cnt = 0;
        rep(i,0,NV) {
            ID[i] = -1;
            vis[i] = -1;
        }
        In[root] = 0;

        rep(i,0,NV) {
            ret += In[i];
            int v = i;
            while(vis[v] != i && ID[v] == -1 && v != root) {
                vis[v] = i;
                v = pre[v];
            }
            if(v != root && ID[v] == -1) {
                for(u = pre[v]; u != v; u = pre[u]) {
                    ID[u] = cnt;
                }
                ID[v] = cnt++;
            }
        }
        if(cnt == 0)    break;
        rep(i,0,NV) {
            if(ID[i] == -1) ID[i] = cnt++;
        }
        rep(i,0,NE) {
            v = E_copy[i].v;
            E_copy[i].u = ID[E_copy[i].u];
            E_copy[i].v = ID[E_copy[i].v];
            if(E_copy[i].u != E_copy[i].v) {
                E_copy[i].cost -= In[v];
            }
        }
        NV = cnt;
        root = ID[root];
    }
    return ret;
}
}

int adj[105][105];
int main() {
	rep(it,0,50000) {
		bumpalloc.reset();
		int n = (rand()%20)+1;
		int density = rand() % 101;
		int r = rand()%n;
		int cnt = 0;
		vector<Edge> edges;
		rep(i,0,n)
			rep(j,0,n){
				if (i==j) continue;
				if (rand() % 100 >= density) continue;
				int weight = rand()%100;
				mit::E[cnt++] = {i,j, weight};
				edges.push_back({i,j,weight});
				adj[i][j] = weight;
			}

		ll ans1 = mit::Directed_MST(r, n, cnt);
		auto pa = dmst(n, r, edges);
		ll ans2 = pa.first;
		assert(ans1 == ans2);

		// Verifying reconstruction:
		if (ans1 != -1) {
			vi par = pa.second;
			if (0) {
				cout << "r = " << r << endl;
				for(auto &x: par) cout << x << ' ';
				cout << endl;
				for(auto &e: edges) {
					cout << e.a << ' ' << e.b << ' ' << e.w << endl;
				}
			}
			ll sum = 0;
			vector<vi> ch(n);
			rep(i,0,n) {
				if (i == r) assert(par[i] == -1);
				else {
					assert(par[i] != -1);
					sum += adj[par[i]][i];
					ch[par[i]].push_back(i);
				}
			}
			assert(sum == ans1);
			vi seen(n), q = {r};
			rep(qi,0,sz(q)) {
				int s = q[qi];
				if (!seen[s]++)
					for(auto &x: ch[s]) q.push_back(x);
			}
			assert(count(all(seen), 0) == 0);
		}
	}
	// Multi-edges, self-loops, negative and huge weights against an
	// exhaustive search over all parent assignments.
	rep(it,0,100000) {
		bumpalloc.reset();
		int n = rand() % 6 + 1, r = rand() % n;
		int m = rand() % 3 ? rand() % (3 * n + 1) : rand() % (n + 1);
		int mode = rand() % 4;
		const ll INF = LLONG_MAX;
		vector<Edge> edges;
		vector<vector<ll>> best(n, vector<ll>(n, INF));
		rep(i,0,m) {
			int a = rand() % n, b = rand() % n;
			if (mode == 3 && rand() % 2) b = a; // extra self-loops
			ll w = mode == 0 ? rand() % 5 : mode == 1 ? rand() % 21 - 10
				: (ll)(rand() % 2001 - 1000) * 1000000000000LL;
			edges.push_back({a, b, w});
			if (a != b) best[a][b] = min(best[a][b], w);
		}
		ll opt = INF;
		vi par(n, -1), cur(n);
		function<void(int, ll)> go = [&](int v, ll sum) {
			if (v == n) {
				rep(i,0,n) {
					int x = i, steps = 0;
					while (x != r && steps++ <= n) x = par[x];
					if (x != r) return;
				}
				opt = min(opt, sum);
				return;
			}
			if (v == r) return go(v + 1, sum);
			rep(p,0,n) if (best[p][v] != INF) {
				par[v] = p;
				go(v + 1, sum + best[p][v]);
			}
		};
		go(0, 0);
		auto pa = dmst(n, r, edges);
		if (opt == INF) {
			assert(pa.first == -1 && pa.second.empty());
			continue;
		}
		assert(pa.first == opt && sz(pa.second) == n);
		ll sum = 0;
		rep(i,0,n) {
			int p = pa.second[i];
			if (i == r) { assert(p == -1); continue; }
			assert(0 <= p && p < n && p != i && best[p][i] != INF);
			sum += best[p][i];
			int x = i, steps = 0;
			while (x != r && steps++ <= n) x = pa.second[x];
			assert(x == r);
		}
		assert(sum == opt);
	}
	cout<<"Tests passed!"<<endl;
	return 0;
}
