#include "../utilities/template.h"

// MinCut.h is a note: after max-flow, the vertices reachable from s in the
// residual graph form the s-side of a min cut. Check that claim with a
// simple Edmonds-Karp against brute force over all s-t cuts.
int main() {
	srand(9);
	rep(it,0,100000) {
		int n = rand() % 8 + 2, m = rand() % 25;
		int s = rand() % n, t = rand() % n;
		if (s == t) t = (s + 1) % n;
		vector<vector<ll>> cap(n, vector<ll>(n));
		rep(i,0,m) { // multi-edges merge, self-loops and zero caps allowed
			int a = rand() % n, b = rand() % n;
			cap[a][b] += rand() % 3 ? rand() % 10 : rand() % 2 * (ll)1e12;
		}
		auto res = cap;
		ll flow = 0;
		vi par(n);
		auto bfs = [&]() {
			fill(all(par), -1); par[s] = s;
			vi q = {s};
			rep(i,0,sz(q)) rep(y,0,n) if (par[y] < 0 && res[q[i]][y] > 0)
				par[y] = q[i], q.push_back(y);
			return par[t] >= 0;
		};
		while (bfs()) {
			ll f = LLONG_MAX;
			for (int x = t; x != s; x = par[x]) f = min(f, res[par[x]][x]);
			for (int x = t; x != s; x = par[x])
				res[par[x]][x] -= f, res[x][par[x]] += f;
			flow += f;
		}
		// par[] now marks the vertices reachable from s
		assert(par[s] >= 0 && par[t] < 0);
		ll cut = 0;
		rep(i,0,n) rep(j,0,n) if (par[i] >= 0 && par[j] < 0) cut += cap[i][j];
		assert(cut == flow);
		ll best = LLONG_MAX;
		rep(mask,0,1 << n) if ((mask >> s & 1) && !(mask >> t & 1)) {
			ll c = 0;
			rep(i,0,n) rep(j,0,n)
				if ((mask >> i & 1) && !(mask >> j & 1)) c += cap[i][j];
			best = min(best, c);
		}
		assert(best == flow);
	}
	cout<<"Tests passed!"<<endl;
}
