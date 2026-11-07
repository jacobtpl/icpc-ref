#include "../utilities/template.h"
#define pb push_back

#include "../../content/graph/FlowDemands.h"
#include "../../content/graph/EdmondsKarp.h"

struct E { int u, v; ll lo, hi; };

// random ll in [lo, hi]
ll rnd(ll lo, ll hi) {
	ll r = ((ll)rand() << 31 ^ rand()) % (hi - lo + 1);
	return lo + r;
}

// Checks the per-edge flows reported by F: bounds, conservation with supplies,
// and net flow `val` out of s / into t (s = t = -1 for a circulation).
void verify(FlowDemands& F, const vector<E>& es, const vector<ll>& sup, int s, int t, ll val) {
	vector<ll> net(sz(sup));
	rep(i,0,sz(es)) {
		ll f = F.flow(i);
		assert(es[i].lo <= f && f <= es[i].hi);
		net[es[i].u] += f; net[es[i].v] -= f;
	}
	rep(v,0,sz(sup)) {
		ll want = sup[v];
		if (v == s) want += val;
		if (v == t) want -= val;
		assert(net[v] == want);
	}
}

// Brute force: enumerate every integer flow within the bounds.
// Returns max value of a feasible s-t flow with value >= 0, or -1 (s = t = -1: 0 if feasible).
ll brute(int n, const vector<E>& es, const vector<ll>& sup, int s, int t) {
	int m = sz(es);
	vector<ll> f(m);
	ll best = -1;
	function<void(int)> go = [&](int i) {
		if (i == m) {
			vector<ll> net(n);
			rep(j,0,m) net[es[j].u] += f[j], net[es[j].v] -= f[j];
			ll val = 0;
			if (s >= 0) {
				val = net[s] - sup[s];
				if (val < 0 || net[t] - sup[t] != -val) return;
			}
			rep(v,0,n) if (v != s && v != t && net[v] != sup[v]) return;
			best = max(best, val);
			return;
		}
		for (f[i] = es[i].lo; f[i] <= es[i].hi; f[i]++) go(i + 1);
	};
	go(0);
	return best;
}

// Independent implementation of the same reduction on top of EdmondsKarp.h.
ll reference(int n, const vector<E>& es, const vector<ll>& sup, int s, int t, ll INF) {
	vector<unordered_map<int, ll>> g(n + 2);
	vector<ll> ex = sup;
	for (auto& e : es) {
		if (e.lo > e.hi) return -1;
		if (e.u != e.v) g[e.u][e.v] += e.hi - e.lo;
		ex[e.v] += e.lo; ex[e.u] -= e.lo;
	}
	if (s >= 0) g[t][s] += INF;
	ll need = 0;
	rep(i,0,n) {
		if (ex[i] > 0) g[n][i] += ex[i], need += ex[i];
		else if (ex[i] < 0) g[i][n + 1] += -ex[i];
	}
	if (edmondsKarp(g, n, n + 1) != need) return -1;
	if (s < 0) return 0;
	return edmondsKarp(g, s, t);
}

bool bounded(int n, const vector<E>& es, int s, int t) {
	vector<vi> adj(n);
	for (auto& e : es) if (e.hi == LLONG_MAX) adj[e.u].pb(e.v);
	vector<bool> vis(n); vis[s] = 1;
	vi q = {s};
	rep(i,0,sz(q)) for (int v : adj[q[i]]) if (!vis[v]) vis[v] = 1, q.pb(v);
	return !vis[t];
}

ll solve(int n, const vector<E>& es, const vector<ll>& sup, int s, int t) {
	FlowDemands F(n);
	for (auto& e : es) assert(F.addEdge(e.u, e.v, e.lo, e.hi) == sz(F.lo) - 1);
	rep(v,0,n) F.ex[v] += sup[v];
	ll val = s < 0 ? (F.circulation() ? 0 : -1) : F.maxFlow(s, t);
	if (val >= 0) verify(F, es, sup, s, t, val);
	return val;
}

void pickST(int n, int& s, int& t) {
	s = rand() % n;
	t = rand() % (n - 1);
	if (t >= s) t++;
}

vector<ll> randSupplies(int n) {
	vector<ll> sup(n);
	if (n >= 2 && rand() % 2) {
		rep(i,0,n) sup[i] = rand() % 5 - 2;
		sup[rand() % n] -= accumulate(all(sup), 0LL);
	}
	return sup;
}

void testTiny() {
	rep(it,0,30000) {
		int n = 1 + rand() % 4, m = rand() % 6;
		vector<E> es;
		rep(i,0,m) {
			ll lo = rand() % 5 - 2;
			ll hi = rand() % 10 == 0 ? lo - 1 : lo + rand() % 4;
			es.pb({rand() % n, rand() % n, lo, hi});
		}
		vector<ll> sup = randSupplies(n);
		int s = -1, t = -1;
		if (n >= 2 && rand() % 4) pickST(n, s, t);
		ll want = brute(n, es, sup, s, t);
		assert(solve(n, es, sup, s, t) == want);
		assert(reference(n, es, sup, s, t, 1000) == want);
	}
}

void testRandom() {
	rep(it,0,20000) {
		int n = 2 + rand() % 9, m = rand() % 26;
		bool big = rand() % 2;
		ll R = big ? (ll)1e12 : 20;
		vector<E> es;
		rep(i,0,m) {
			ll lo = rnd(-R, R);
			ll hi = rand() % 15 == 0 ? lo - 1 - rnd(0, R) : lo + rnd(0, R);
			es.pb({rand() % n, rand() % n, lo, hi});
		}
		vector<ll> sup = randSupplies(n);
		if (big) for (ll& x : sup) x *= (ll)1e11;
		int s = -1, t = -1;
		if (rand() % 4) pickST(n, s, t);
		ll want = reference(n, es, sup, s, t, (ll)1e17);
		assert(solve(n, es, sup, s, t) == want);
	}
}

// hi = LLONG_MAX edges (lo >= 0), compared against the same graph with hi capped.
void testInfinite() {
	int done = 0;
	while (done < 5000) {
		int n = 2 + rand() % 7, m = rand() % 16;
		vector<E> es, capped;
		rep(i,0,m) {
			ll lo = rnd(-(ll)1e12, (ll)1e12);
			ll hi = lo + rnd(0, (ll)1e12);
			if (rand() % 4 == 0) lo = rnd(0, (ll)1e12), hi = LLONG_MAX;
			es.pb({rand() % n, rand() % n, lo, hi});
			capped.pb({es[i].u, es[i].v, lo, hi == LLONG_MAX ? (ll)1e17 : hi});
		}
		vector<ll> sup(n);
		int s = -1, t = -1;
		if (rand() % 4) {
			pickST(n, s, t);
			if (!bounded(n, es, s, t)) continue;
		}
		done++;
		ll want = reference(n, capped, sup, s, t, (ll)1e17);
		assert(solve(n, capped, sup, s, t) == want);
		assert(solve(n, es, sup, s, t) == want);
	}
}

void testEdgeCases() {
	// no nodes / no edges
	assert(solve(0, {}, {}, -1, -1) == 0);
	assert(solve(1, {}, {0}, -1, -1) == 0);
	assert(solve(2, {}, {0, 0}, 0, 1) == 0);
	// n = 1: only self-loops
	assert(solve(1, {{0, 0, 3, 5}}, {0}, -1, -1) == 0);
	assert(solve(1, {{0, 0, -4, -4}}, {0}, -1, -1) == 0);
	assert(solve(1, {{0, 0, 5, 3}}, {0}, -1, -1) == -1);
	// lo > hi anywhere is infeasible
	assert(solve(2, {{0, 1, 5, 3}, {1, 0, 0, 10}}, {0, 0}, -1, -1) == -1);
	assert(solve(2, {{0, 1, 5, 3}}, {0, 0}, 0, 1) == -1);
	assert(solve(3, {{0, 1, 0, 10}, {1, 2, 0, 10}, {1, 1, 1, 0}}, {0, 0, 0}, 0, 2) == -1);
	// unsatisfiable lower bound (no way back)
	assert(solve(2, {{0, 1, 1, 5}}, {0, 0}, -1, -1) == -1);
	assert(solve(2, {{0, 1, 1, 5}, {1, 0, 0, 3}}, {0, 0}, -1, -1) == 0);
	// forced flow 1e18 along a chain, LLONG_MAX return edge
	assert(solve(3, {{0, 1, (ll)1e18, (ll)1e18}, {1, 2, (ll)1e18, (ll)1e18}}, {0, 0, 0}, 0, 2) == (ll)1e18);
	assert(solve(3, {{0, 1, (ll)1e18, (ll)1e18}, {1, 2, (ll)1e18, (ll)1e18}}, {0, 0, 0}, -1, -1) == -1);
	assert(solve(3, {{0, 1, (ll)1e18, (ll)1e18}, {1, 2, (ll)1e18, (ll)1e18}, {2, 0, 0, LLONG_MAX}}, {0, 0, 0}, -1, -1) == 0);
	// a forced s-t flow larger than what can return is fine for maxFlow, infeasible as a circulation
	assert(solve(2, {{0, 1, 7, 7}, {1, 0, 0, 3}}, {0, 0}, 0, 1) == 7);
	assert(solve(2, {{0, 1, 7, 7}, {1, 0, 0, 3}}, {0, 0}, -1, -1) == -1);
	// negative lower bound: flow may go backwards
	assert(solve(2, {{0, 1, -5, 5}}, {0, 0}, 1, 0) == 5);
	assert(solve(2, {{0, 1, -5, -2}}, {0, 0}, 0, 1) == -1);
	// supplies: 0 must push 4 out, 2 must absorb 4, through 1
	assert(solve(3, {{0, 1, 0, 10}, {1, 2, 0, 3}}, {4, 0, -4}, -1, -1) == -1);
	assert(solve(3, {{0, 1, 0, 10}, {1, 2, 0, 4}}, {4, 0, -4}, -1, -1) == 0);
	// supply on s is not counted in the value
	assert(solve(2, {{0, 1, 0, 10}}, {4, -4}, 0, 1) == 6);
	assert(solve(2, {{0, 1, 0, 10}}, {-4, 4}, 0, 1) == 14);
	// parallel edges and both directions, large values
	assert(solve(2, {{0, 1, 0, (ll)4e18}, {0, 1, 0, (ll)4e18}, {1, 0, 0, (ll)1e18}}, {0, 0}, 0, 1) == (ll)8e18);
	// edge index is returned and flow() is per edge
	{
		FlowDemands F(3);
		int a = F.addEdge(0, 1, 2, 2), b = F.addEdge(1, 2, 0, 5), c = F.addEdge(0, 2, 0, 1);
		assert(a == 0 && b == 1 && c == 2);
		assert(F.maxFlow(0, 2) == 3);
		assert(F.flow(a) == 2 && F.flow(b) == 2 && F.flow(c) == 1);
	}
}

int main() {
	testEdgeCases();
	testTiny();
	testRandom();
	testInfinite();
	cout << "Tests passed!" << endl;
}
