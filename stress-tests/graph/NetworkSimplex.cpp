#include "../utilities/template.h"
#include <unistd.h>

// NetworkSimplex.h relies on these helpers from content/contest/template.cpp
#define mp make_pair
#define pb push_back
template<class A, class B> bool ckmin(A& a, const B& b) { return b < a ? a = b, 1 : 0; }

#include "../../content/graph/NetworkSimplex.h"

struct Ed { int a, b; ll l, u, c; };
const ll NONE = LLONG_MIN;

// Oracle 1: enumerate every integral flow vector.
ll brute(int n, const vector<Ed>& es, const vector<ll>& sup) {
	int m = sz(es);
	vector<ll> f(m), bal(n);
	ll best = NONE;
	function<void(int, ll)> go = [&](int i, ll cost) {
		if (i == m) {
			rep(v,0,n) if (bal[v] != sup[v]) return;
			if (best == NONE || cost < best) best = cost;
			return;
		}
		const Ed& e = es[i];
		for (ll x = e.l; x <= e.u; x++) {
			bal[e.a] += x, bal[e.b] -= x;
			go(i + 1, cost + x * e.c);
			bal[e.a] -= x, bal[e.b] += x;
		}
	};
	go(0, 0);
	return best;
}

// Oracle 2: saturate negative edges, then successive shortest paths (Bellman-Ford).
ll ssp(int n, const vector<Ed>& es, const vector<ll>& sup) {
	struct E { int to; ll cap, cost; };
	vector<E> E_; vector<vi> g(n + 2);
	auto add = [&](int a, int b, ll cap, ll cost) {
		g[a].push_back(sz(E_)); E_.push_back({b, cap, cost});
		g[b].push_back(sz(E_)); E_.push_back({a, 0, -cost});
	};
	vector<ll> ex(sup);
	ll cost = 0, need = 0;
	for (auto& e : es) {
		if (e.c < 0) add(e.b, e.a, e.u - e.l, -e.c), ex[e.b] += e.u, ex[e.a] -= e.u, cost += e.u * e.c;
		else add(e.a, e.b, e.u - e.l, e.c), ex[e.b] += e.l, ex[e.a] -= e.l, cost += e.l * e.c;
	}
	int S = n, T = n + 1;
	rep(v,0,n) {
		if (ex[v] > 0) add(S, v, ex[v], 0), need += ex[v];
		if (ex[v] < 0) add(v, T, -ex[v], 0);
	}
	ll sum = 0; rep(v,0,n) sum += ex[v];
	if (sum != 0) return NONE;
	for (;;) {
		vector<ll> d(n + 2, LLONG_MAX); vi pre(n + 2, -1);
		d[S] = 0;
		for (bool ch = 1; ch;) {
			ch = 0;
			rep(v,0,n+2) if (d[v] != LLONG_MAX) for (int id : g[v])
				if (E_[id].cap > 0 && d[v] + E_[id].cost < d[E_[id].to])
					d[E_[id].to] = d[v] + E_[id].cost, pre[E_[id].to] = id, ch = 1;
		}
		if (d[T] == LLONG_MAX) break;
		ll f = LLONG_MAX;
		for (int v = T; v != S; v = E_[pre[v] ^ 1].to) f = min(f, E_[pre[v]].cap);
		for (int v = T; v != S; v = E_[pre[v] ^ 1].to) E_[pre[v]].cap -= f, E_[pre[v] ^ 1].cap += f;
		cost += f * d[T], need -= f;
	}
	return need ? NONE : cost;
}

// Runs NetworkSimplex and validates the flow it leaves in E. Returns NONE if infeasible.
ll run(int n, const vector<Ed>& es, const vector<ll>& sup) {
	NetworkSimplex ns; ns.init(n);
	for (auto& e : es) ns.ae(e.a, e.b, e.l, e.u, e.c);
	rep(v,0,n) ns.B[v] += sup[v];
	i128 r;
	try { r = ns.solve(); } catch (int) { return NONE; }
	vector<ll> bal(n); i128 cost = 0;
	rep(i,0,sz(es)) {
		ll f = ns.E[2 * i].flow + es[i].l;
		assert(es[i].l <= f && f <= es[i].u);
		assert(ns.E[2 * i + 1].flow == -ns.E[2 * i].flow);
		bal[es[i].a] += f, bal[es[i].b] -= f;
		cost += (i128)f * es[i].c;
	}
	rep(v,0,n) assert(bal[v] == sup[v]);
	assert(cost == r);
	return (ll)r;
}

mt19937 rng(12345);
int rnd(int lo, int hi) { return lo + (int)(rng() % (unsigned)(hi - lo + 1)); }

// kind 0: circulation only; kind 1: random supplies summing to 0; kind 2: arbitrary supplies
void gen(int n, int m, int capMax, int costMax, bool lower, bool loops, int kind,
		vector<Ed>& es, vector<ll>& sup) {
	es.clear(); sup.assign(n, 0);
	rep(i,0,m) {
		int a = rnd(0, n - 1), b = rnd(0, n - 1);
		if (!loops && a == b) { if (n == 1) continue; b = (a + 1) % n; }
		int l = lower ? rnd(-capMax, capMax) : 0, u = l + rnd(0, capMax);
		es.push_back({a, b, l, u, rnd(-costMax, costMax)});
	}
	if (kind == 1) rep(i,0,rnd(0, 3)) {
		int a = rnd(0, n - 1), b = rnd(0, n - 1), x = rnd(0, capMax);
		sup[a] += x, sup[b] -= x;
	}
	if (kind == 2) rep(v,0,n) sup[v] = rnd(-1, 1);
}

int main() {
	vector<Ed> es; vector<ll> sup;
	// empty / trivial instances
	assert(run(1, {}, {0}) == 0);
	assert(run(1, {}, {1}) == NONE);
	assert(run(2, {}, {1, -1}) == NONE);
	assert(run(2, {{0, 1, 0, 5, 3}}, {2, -2}) == 6);
	assert(run(2, {{0, 1, 2, 5, 3}, {1, 0, 0, 5, -1}}, {0, 0}) == 4);
	assert(run(2, {{0, 1, 0, 5, -3}, {1, 0, 0, 4, 1}}, {0, 0}) == -8);
	assert(run(1, {{0, 0, 0, 3, -2}}, {0}) == -6); // self-loops
	assert(run(1, {{0, 0, 1, 3, 2}}, {0}) == 2);
	assert(run(2, {{0, 1, 0, 0, -5}, {1, 0, 0, 7, 1}}, {0, 0}) == 0); // zero capacity
#ifdef NS_EMPTY_GRAPH
	// n = 0 reads E[0] of an empty vector (segfault); opt-in, see audit notes
	assert(run(0, {}, {}) == 0);
#endif
	bool loops = 1;
	// exhaustive oracle on tiny instances
	rep(it,0,150000) {
		int n = rnd(1, 4), m = rnd(0, 5);
		gen(n, m, rnd(1, 2), rnd(0, 3), it % 2, loops && it % 3 == 0, it % 5 % 3, es, sup);
		ll a = run(n, es, sup), b = brute(n, es, sup), c = ssp(n, es, sup);
		assert(a == b); assert(b == c);
	}
	// SSP oracle on larger instances
	rep(it,0,6000) {
		int n = rnd(1, 12), m = rnd(0, 40);
		gen(n, m, rnd(1, 20), it % 4 ? rnd(0, 50) : 1, it % 2, loops && it % 3 == 0, it % 5 % 3, es, sup);
		assert(run(n, es, sup) == ssp(n, es, sup));
	}
	rep(it,0,150) {
		int n = rnd(20, 60), m = rnd(n, 6 * n);
		gen(n, m, 1000, 1000, it % 2, 0, it % 2, es, sup);
		assert(run(n, es, sup) == ssp(n, es, sup));
	}
	// highly degenerate instances (many ties / zero pivots): must not cycle
	rep(it,0,300) {
		int n = rnd(10, 50), m = rnd(n, 8 * n);
		gen(n, m, 1, it % 3 ? 1 : 0, 0, 1, it % 2, es, sup);
		assert(run(n, es, sup) == ssp(n, es, sup));
	}
	// s-t min-cost max-flow through a very negative t->s edge: mostly degenerate pivots.
	// With an arbitrary leaving-arc rule this stalls (> 3e6 pivots, minutes); it needs ~2e3.
	{
		mt19937 r2(1);
		int n = 500; ll BIG = (ll)2e9, tot = 0;
		es.clear(); sup.assign(n, 0);
		rep(i,0,10 * n) {
			int a = (int)(r2() % n), b = (int)(r2() % n);
			if (a == b) continue;
			es.push_back({a, b, 0, (int)(r2() % 10) + 1, (int)(r2() % 1000)});
		}
		for (auto& e : es) if (e.a == 0) tot += e.u;
		es.push_back({n - 1, 0, 0, tot, -BIG});
		alarm(20); // fail instead of hanging
		ll a = run(n, es, sup);
		alarm(0);
		assert(a == ssp(n, es, sup));
	}
	// lower bound * cost exceeds 64 bits; the result is an i128
	{
		ll L = (ll)4e9, C = (ll)4e9;
		NetworkSimplex ns; ns.init(2);
		ns.ae(0, 1, L, L, C); ns.ae(1, 0, 0, L, 0);
		assert(ns.solve() == (i128)L * C);
	}
	// large costs and capacities: answer needs more than 64 bits to be safe
	{
		ll C = (ll)1e12, U = (ll)1e12;
		NetworkSimplex ns; ns.init(3);
		ns.ae(0, 1, 0, U, -C); ns.ae(1, 2, 0, U, -C); ns.ae(2, 0, 0, U, -C);
		assert(ns.solve() == -(i128)3 * C * U);
	}
	cout << "Tests passed!" << endl;
}
