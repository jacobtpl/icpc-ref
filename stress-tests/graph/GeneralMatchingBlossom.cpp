#include "../utilities/template.h"
#define pb push_back // from content/contest/template.cpp

#include "../../content/graph/GeneralMatchingBlossom.h"

mt19937 rng(2024);
int ri(int a, int b) { return (int)(rng() % (unsigned)(b - a + 1)) + a; }

// Maximum matching size by bitmask DP, ignoring vertex `skip`.
int brute(int n, const vector<vi>& g, int skip) {
	vi dp(1 << n);
	rep(m,1,1<<n) {
		int i = __builtin_ctz(m);
		dp[m] = dp[m ^ (1 << i)];
		if (i != skip) rep(j,i+1,n) if ((m >> j & 1) && g[i][j] && j != skip)
			dp[m] = max(dp[m], dp[m ^ (1 << i) ^ (1 << j)] + 1);
	}
	return dp[(1 << n) - 1];
}

// Independent O(N^3)-ish oracle for larger graphs: simple augmenting path
// blossom algorithm (classic contraction-based Edmonds with base[] arrays).
struct Edmonds {
	int n; vector<vi> g; vi mt, p, base; vector<char> used, bl;
	Edmonds(int n) : n(n), g(n), mt(n, -1), p(n), base(n) {}
	int lca(int a, int b) {
		vector<char> seen(n);
		for (;;) { a = base[a]; seen[a] = 1; if (mt[a] == -1) break; a = p[mt[a]]; }
		for (;;) { b = base[b]; if (seen[b]) return b; b = p[mt[b]]; }
	}
	void mark(int v, int b, int ch) {
		while (base[v] != b) {
			bl[base[v]] = bl[base[mt[v]]] = 1;
			p[v] = ch; ch = mt[v]; v = p[mt[v]];
		}
	}
	int path(int root) {
		used.assign(n, 0); p.assign(n, -1);
		rep(i,0,n) base[i] = i;
		used[root] = 1; queue<int> q; q.push(root);
		while (!q.empty()) {
			int v = q.front(); q.pop();
			for (int to : g[v]) {
				if (base[v] == base[to] || mt[v] == to) continue;
				if (to == root || (mt[to] != -1 && p[mt[to]] != -1)) {
					int cb = lca(v, to); bl.assign(n, 0);
					mark(v, cb, to); mark(to, cb, v);
					rep(i,0,n) if (bl[base[i]]) {
						base[i] = cb;
						if (!used[i]) used[i] = 1, q.push(i);
					}
				} else if (p[to] == -1) {
					p[to] = v;
					if (mt[to] == -1) return to;
					used[mt[to]] = 1; q.push(mt[to]);
				}
			}
		}
		return -1;
	}
	int solve() {
		int res = 0;
		rep(i,0,n) if (mt[i] == -1) {
			int v = path(i);
			if (v != -1) res++;
			while (v != -1) { int pv = p[v], ppv = mt[pv]; mt[v] = pv; mt[pv] = v; v = ppv; }
		}
		return res;
	}
};

// Runs MaxMatching and validates the returned matching; returns its size.
int run(int n, const vector<pii>& ed, MaxMatching& M) {
	M.init(n);
	set<pii> es;
	for (pii e : ed) {
		M.ae(e.first + 1, e.second + 1);
		es.insert(e), es.insert({e.second, e.first});
	}
	int ans = M.solve(), cnt = 0;
	assert(M.mate[0] == 0);
	rep(i,1,n+1) if (int x = M.mate[i]) {
		assert(1 <= x && x <= n && x != i && M.mate[x] == i);
		assert(es.count({i - 1, x - 1}));
		cnt++;
	}
	assert(cnt == 2 * ans);
	return ans;
}

void testSmall(int maxn, int iters, bool loops) {
	rep(it,0,iters) {
		int n = ri(0, maxn), m = n ? ri(0, n * n) : 0, dens = ri(0, 100);
		vector<vi> g(n, vi(n)); vector<pii> ed;
		rep(e,0,m) { // multi-edges, optionally self-loops
			int a = ri(0, n-1), b = ri(0, n-1);
			if (ri(0, 99) >= dens || (a == b && !loops)) continue;
			ed.push_back({a, b});
			if (a != b) g[a][b] = g[b][a] = 1;
		}
		MaxMatching M;
		int ans = run(n, ed, M), best = brute(n, g, -1);
		assert(ans == best);
		// documented: white[v] == 0  =>  v is in every maximum matching
		// (the converse also holds in practice, so check both directions)
		rep(v,0,n)
			assert(!M.white[v + 1] == (brute(n, g, v) < best));
	}
}

void testLarge(int maxn, int iters) {
	rep(it,0,iters) {
		int n = ri(1, maxn), type = ri(0, 3), m;
		if (type == 0) m = ri(0, 2 * n);          // sparse: many blossoms
		else if (type == 1) m = ri(0, n * n / 4); // dense
		else m = n;
		vector<pii> ed;
		rep(e,0,m) {
			int a = ri(0, n-1), b = ri(0, n-1);
			if (type == 2) a = e, b = (e + 1) % n;        // one big cycle
			if (type == 3) a = e, b = e % 3 ? e - 1 : max(e - 3, 0); // triangles chain
			if (type == 3 && e % 3 == 2) ed.push_back({e, e - 2});
			if (a != b) ed.push_back({a, b});
		}
		shuffle(all(ed), rng);
		MaxMatching M;
		Edmonds E(n);
		for (pii e : ed) E.g[e.first].push_back(e.second), E.g[e.second].push_back(e.first);
		assert(run(n, ed, M) == E.solve());
	}
}

int main() {
	testSmall(1, 100, 1);
	testSmall(5, 50000, 0);
	testSmall(5, 20000, 1);
	testSmall(9, 30000, 0);
	testSmall(9, 10000, 1);
	testSmall(13, 1000, 0);
	testLarge(30, 20000);
	testLarge(300, 300);
	cout<<"Tests passed!"<<endl;
}
