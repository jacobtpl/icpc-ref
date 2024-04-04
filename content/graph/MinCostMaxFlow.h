/**
 * Author: Siyong, modified from Kactl Implementation
 * Date: April 4, 2024
 * Description: Min-cost max-flow. All capacities are 0. Flows are initialized to be negative.
 * Status: Tested on kattis:mincostmaxflow
 * Time: Originally O(E^2)
 */

// #include <bits/extc++.h>
struct MCMF {
	typedef int C; typedef int F; typedef ll R;
	C const INFC = numeric_limits<C>::max() / 4;
	F const INFF = numeric_limits<F>::max() / 4;

	struct Edge {int n; F flow; C cost; size_t rev;};
	int N;
	vector<vector<Edge>> ed;
	vector<C> dist, pi;
	vector<F> amt;
	vector<size_t> par;

	MCMF(int N) :
		N(N), ed(N), dist(N), pi(N), par(N), amt(N) {}

	void addEdge(int u, int v, F f, C c) {
		ed[u].emplace_back(v, -f, c, ed[v].size());
		ed[v].emplace_back(u, 0, -c, ed[u].size()-1);
	}

	void path(int s) { // 417ab0
		fill(all(amt), 0); amt[s] = INFF;
		fill(all(dist), INFC); dist[s] = 0;

		__gnu_pbds::priority_queue<pair<C, int>> q;
		vector<decltype(q)::point_iterator> its(N);
		q.push({0, s});

		while (!q.empty()) {
			s = q.top().second; q.pop();
			C di = dist[s] + pi[s];
			for (auto [n, f, c, rev]: ed[s])
				if (f < 0 && ckmin(dist[n], di + c - pi[n])) {
					par[n] = rev; amt[n] = min(amt[s], -f);
					if(its[n] == q.end()) its[n] = q.push({-dist[n], n});
					else q.modify(its[n], {-dist[n], n});
				}
		}
		rep(i,0,N) pi[i] = amt[i] ? pi[i] + dist[i] : INFC;
	}

	pair<F, R> maxflow(int s, int t) { // 2126d0
		F totflow = 0; R totcost = 0;
		while (path(s), amt[t]) {
			F fl = amt[t]; totflow += fl;
			for(int n = t;n != s;) {
				auto &[p, f, c, rev] = ed[n][par[n]];
				f -= fl; ed[p][rev].flow += fl; n = p;
				totcost -= (R) fl * c;//OR +=(R)fl*ed[p][rev].cost
			}
		}
		return {totflow, totcost};
	}

	// If some costs can be negative, call this before maxflow:
	void setpi(int s) { // (otherwise, leave this out)
		fill(all(pi), INFC); pi[s] = 0;
		int it = N, ch = 1;
		while (ch-- && it--)
			rep(i,0,N) if (pi[i] != INFC)
				for (auto [to, f, c, _]: ed[i])
					if (f < 0 && ckmin(pi[to], pi[i] + c)) ch = 1;
		assert(it >= 0); // negative cost cycle
	}
}; // 0b54fa without setpi; 88d7c5 with setpi
