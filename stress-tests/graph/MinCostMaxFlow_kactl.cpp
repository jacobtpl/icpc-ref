#include "../utilities/template.h"

#include <bits/extc++.h>
#define setpi dummy(){} bool setpi
#undef assert
#define assert(x) return x
#include "../../content/graph/MinCostMaxFlow_kactl.h"
#undef assert
#undef setpi
#include <cassert>

// Independent oracle: successive shortest paths with Bellman-Ford.
struct RefMCMF {
	struct E { int u, v; ll cap, cost; };
	int n; vector<E> es;
	RefMCMF(int n) : n(n) {}
	void addEdge(int u, int v, ll cap, ll cost) {
		if (u == v) return;
		es.push_back({u, v, cap, cost});
		es.push_back({v, u, 0, -cost});
	}
	bool negCycle(int s) { // negative cycle reachable from s?
		vector<ll> d(n, LLONG_MAX / 4); d[s] = 0;
		rep(it,0,n+1) {
			bool ch = 0;
			for (E& e : es) if (e.cap > 0 && d[e.u] < LLONG_MAX / 4
					&& d[e.u] + e.cost < d[e.v]) d[e.v] = d[e.u] + e.cost, ch = 1;
			if (!ch) return false;
		}
		return true;
	}
	pair<ll, ll> maxflow(int s, int t) {
		ll fl = 0, co = 0;
		const ll inf = LLONG_MAX / 4;
		for (;;) {
			vector<ll> d(n, inf); vi pe(n, -1);
			d[s] = 0;
			rep(it,0,n) rep(i,0,sz(es)) {
				E& e = es[i];
				if (e.cap > 0 && d[e.u] < inf && d[e.u] + e.cost < d[e.v])
					d[e.v] = d[e.u] + e.cost, pe[e.v] = i;
			}
			if (d[t] == inf) break;
			ll f = inf;
			for (int x = t; x != s; x = es[pe[x]].u) f = min(f, es[pe[x]].cap);
			for (int x = t; x != s; x = es[pe[x]].u)
				es[pe[x]].cap -= f, es[pe[x] ^ 1].cap += f;
			fl += f, co += f * d[t];
		}
		return {fl, co};
	}
};

// mode 0: costs >= 0; mode 1: negative costs, no negative cycles;
// mode 2: arbitrary small costs, instances where setpi reports a negative
// cycle are checked against the oracle and skipped.
// No double edges (documented), but antiparallel edges, self-loops and
// zero capacities are generated.
ll done[3];
void testRef(int its, int maxN, int maxM, ll maxCap, int maxCost, int mode) {
	rep(it,0,its) {
		int N = rand() % maxN + 1, M = rand() % (maxM + 1);
		int S = rand() % N, T = rand() % N;
		if (S == T) { if (N == 1) continue; T = (S + 1) % N; }
		MCMF mcmf(N); RefMCMF ref(N);
		vi pot(N);
		if (mode == 1) rep(i,0,N) pot[i] = rand() % (maxCost + 1);
		struct Ad { int u, v; ll cap, cost; };
		vector<Ad> ads;
		set<pii> used;
		rep(i,0,M) {
			int u = rand() % N, v = rand() % N;
			if (mode && u == v) continue;
			if (!used.insert({u, v}).second) continue;
			ll cap = rand() % 4 == 0 ? 0 :
				(ll)((unsigned long long)rand() * rand() % (unsigned long long)maxCap) + 1;
			if (rand() % 8 == 0) cap = maxCap;
			ll cost = rand() % (maxCost + 1) + pot[u] - pot[v];
			if (mode == 2) cost = rand() % 11 - 3;
			ads.push_back({u, v, cap, cost});
			mcmf.addEdge(u, v, cap, cost);
			ref.addEdge(u, v, cap, cost);
		}
		if (mode || rand() % 4 == 0) {
			bool ok = mcmf.setpi(S);
			assert(ok == !ref.negCycle(S));
			if (!ok) continue;
		}
		auto pa = mcmf.maxflow(S, T);
		auto pb = ref.maxflow(S, T);
		assert(pa == pb);
		// "To obtain the actual flow, look at positive values only."
		vector<ll> bal(N); ll co = 0;
		for (auto& a : ads) if (a.u != a.v) {
			ll f = mcmf.flow[a.u][a.v];
			assert(0 <= f && f <= a.cap);
			bal[a.u] -= f, bal[a.v] += f, co += f * a.cost;
		}
		rep(i,0,N) assert(bal[i] == (i == S ? -pa.first : i == T ? pa.first : 0));
		assert(co == pa.second);
		done[mode]++;
	}
}

int main() {
	srand(4);
	testRef(200000, 6, 12, 5, 6, 0);
	testRef(200000, 6, 12, 5, 6, 1);
	testRef(300000, 7, 16, 50, 0, 2);
	testRef(20000, 12, 80, (ll)1e12, 1000000, 0);
	testRef(20000, 12, 80, (ll)1e12, 1000000, 1);
	rep(i,0,3) assert(done[i] > 10000);
	cout<<"Tests passed!"<<endl;
}
