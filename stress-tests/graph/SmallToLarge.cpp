#include "../utilities/template.h"
#include "../utilities/genTree.h"

#include "../../content/graph/SmallToLarge.h"

mt19937 rng(1234);

// Brute force: subtree vertex sets of every vertex reachable from r.
void bruteSub(vector<vi>& g, int v, int p, vector<set<int>>& sub, vi& reach) {
	reach[v] = 1;
	sub[v].insert(v);
	for (int y : g[v]) if (y != p) {
		bruteSub(g, y, v, sub, reach);
		sub[v].insert(all(sub[y]));
	}
}

// Checks every guarantee the header documents on the given tree/root.
void check(vector<vi>& g, int r, const vi& col, int C) {
	int n = sz(g);
	vector<set<int>> sub(n);
	vi reach(n);
	bruteSub(g, r, -1, sub, reach);
	int nReach = 0;
	rep(v,0,n) nReach += reach[v];
	int maxAdds = 1;
	while ((1 << maxAdds) <= nReach) maxAdds++; // 1 + floor(log2 N)

	set<int> active;
	vi adds(n), answered(n), cnt(C), res(n, -1);
	int d = 0;
	sack(g, r, [&](int v) {
		assert(reach[v]);
		assert(active.insert(v).second);
		assert(++adds[v] <= maxAdds);
		d += !cnt[col[v]]++;
	}, [&](int v) {
		assert(active.erase(v));
		d -= !--cnt[col[v]];
	}, [&](int v) {
		assert(reach[v]);
		assert(!answered[v]++);
		assert(active == sub[v]);
		res[v] = d;
	});
	assert(active.empty() && d == 0);
	rep(v,0,n) {
		assert(answered[v] == reach[v]);
		if (reach[v]) {
			set<int> cs;
			for (int u : sub[v]) cs.insert(col[u]);
			assert(res[v] == sz(cs));
		}
	}
}

vector<vi> build(int n, const vector<pii>& edges, bool directed, int r) {
	vector<vi> g(n);
	if (directed) { // orient edges away from r
		vector<vi> u(n);
		for (auto e : edges) u[e.first].push_back(e.second), u[e.second].push_back(e.first);
		vi q{r}, seen(n); seen[r] = 1;
		rep(i,0,sz(q)) for (int y : u[q[i]]) if (!seen[y])
			seen[y] = 1, g[q[i]].push_back(y), q.push_back(y);
	} else {
		for (auto e : edges) g[e.first].push_back(e.second), g[e.second].push_back(e.first);
	}
	for (auto& a : g) shuffle(all(a), rng);
	return g;
}

vector<pii> genShape(int n, int shape) {
	vector<pii> e;
	if (shape == 0) return genRandomTree(n);
	rep(i,1,n) {
		int p = shape == 1 ? i - 1 : shape == 2 ? 0 : shape == 3 ? (i - 1) / 2
			: shape == 4 ? max(0, i - 1 - rand() % 3) : rand() % i;
		e.emplace_back(p, i);
	}
	return e;
}

void randomCase(int n, int extra) {
	int shape = rand() % 6;
	auto edges = n >= 2 ? genShape(n, shape) : vector<pii>();
	vi perm(n); // random relabelling
	iota(all(perm), 0);
	shuffle(all(perm), rng);
	for (auto& e : edges) e = {perm[e.first], perm[e.second]};
	int r = rand() % n;
	bool directed = rand() % 2;
	int N = n + extra; // extra vertices unreachable from r
	auto g = build(N, edges, directed, r);
	rep(i,n,N) if (i > n && rand() % 2) { // link extras among themselves
		int j = n + rand() % (i - n);
		g[j].push_back(i);
		if (!directed) g[i].push_back(j);
	}
	int C = 1 + rand() % 5;
	vi col(N);
	for (int& c : col) c = rand() % C;
	check(g, r, col, C);
}

int main() {
	// n = 1 edge cases
	{
		vector<vi> g(1);
		vi col{0};
		check(g, 0, col, 1);
	}
	{ // n = 1 reachable, one isolated unreachable vertex
		vector<vi> g(2);
		vi col{0, 0};
		check(g, 1, col, 1);
	}
	// n = 2 directed / undirected, both roots
	rep(dir,0,2) rep(r,0,2) {
		auto g = build(2, {{0, 1}}, dir, r);
		vi col{0, 1};
		check(g, r, col, 2);
	}
	rep(it,0,20000) randomCase(1 + rand() % 12, rand() % 4);
	rep(it,0,2000) randomCase(1 + rand() % 60, rand() % 4);
	rep(it,0,100) randomCase(1 + rand() % 300, rand() % 2);

	// Large case: O(N log N) adds, correct distinct-colour counts vs. small-to-large set merge.
	{
		int n = 200000, C = 1000;
		auto edges = genRandomTree(n);
		auto g = build(n, edges, false, 0);
		vi col(n);
		for (int& c : col) c = rand() % C;
		vi cnt(C), res(n), order, par(n, -1), q{0};
		rep(i,0,sz(q)) for (int y : g[q[i]]) if (y != par[q[i]])
			par[y] = q[i], q.push_back(y);
		vector<set<int>> st(n);
		vi want(n);
		for (int i = n; i--;) {
			int v = q[i];
			st[v].insert(col[v]);
			want[v] = sz(st[v]);
			if (par[v] >= 0) {
				int p = par[v];
				if (sz(st[p]) < sz(st[v])) swap(st[p], st[v]);
				st[p].insert(all(st[v]));
				set<int>().swap(st[v]);
			}
		}
		ll adds = 0;
		int d = 0;
		sack(g, 0, [&](int v) { adds++; d += !cnt[col[v]]++; },
			[&](int v) { d -= !--cnt[col[v]]; }, [&](int v) { res[v] = d; });
		assert(res == want);
		assert(adds <= (ll)n * 19);
		assert(d == 0);
	}
	cout << "Tests passed!" << endl;
}
