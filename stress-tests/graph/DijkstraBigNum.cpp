// Dijkstra with a non-numeric distance type: CF 464E "The Classic Problem".
// Edge weights are 2^x (x up to 1e5); a distance is the root of a persistent
// segment tree over bit positions, compared by hash. Checked against ll brute force.
#include "../utilities/template.h"

#include "../../content/graph/Dijkstra.h"

typedef unsigned long long ull;
const ull HMOD = (1ULL << 61) - 1;
ull mulmod(ull a, ull b) {
	__uint128_t c = (__uint128_t)a * b;
	ull r = (ull)(c >> 61) + (ull)(c & HMOD);
	return r >= HMOD ? r - HMOD : r;
}
const int B = 100000 + 20; // bit positions [0, B)
struct BigBits {
	struct Nd { int l, r, cnt; ull h; };
	vector<Nd> t;
	vector<ull> pw;
	int one = -1;
	BigBits() : t(1, {0, 0, 0, 0}), pw(B + 1, 1) {
		ull base = 1000003 + (ull)rand() % 1000000;
		rep(i,1,B+1) pw[i] = mulmod(pw[i-1], base);
	}
	int mk(int l, int r, int len) {
		int lenL = len - len / 2;
		t.push_back({l, r, t[l].cnt + t[r].cnt,
			(t[l].h + mulmod(t[r].h, pw[lenL])) % HMOD});
		return sz(t) - 1;
	}
	int firstZero(int v, int lo, int hi, int x) {
		if (hi <= x || t[v].cnt == hi - lo) return -1;
		if (hi - lo == 1) return lo;
		int mid = lo + (hi - lo + 1) / 2;
		int a = firstZero(t[v].l, lo, mid, x);
		return a != -1 ? a : firstZero(t[v].r, mid, hi, x);
	}
	int upd(int v, int lo, int hi, int x, int p) { // bits [x,p) := 0, bit p := 1
		if (hi <= x || lo > p) return v;
		if (hi - lo == 1) {
			if (lo != p) return 0;
			if (one < 0) { t.push_back({0, 0, 1, 1}); one = sz(t) - 1; }
			return one;
		}
		if (x <= lo && hi <= p) return 0;
		int mid = lo + (hi - lo + 1) / 2;
		return mk(upd(t[v].l, lo, mid, x, p), upd(t[v].r, mid, hi, x, p), hi - lo);
	}
	int add(int v, int x) { return upd(v, 0, B, x, firstZero(v, 0, B, x)); }
	bool less(int a, int b) {
		int lo = 0, hi = B;
		while (hi - lo > 1) {
			if (t[a].h == t[b].h) return false;
			int mid = lo + (hi - lo + 1) / 2;
			if (t[t[a].r].h != t[t[b].r].h) a = t[a].r, b = t[b].r, lo = mid;
			else a = t[a].l, b = t[b].l, hi = mid;
		}
		return t[a].cnt < t[b].cnt;
	}
	ll value(int v, int lo, int hi) { // exact value, only for small tests
		if (!t[v].cnt) return 0;
		if (hi - lo == 1) return 1LL << lo;
		int mid = lo + (hi - lo + 1) / 2;
		return value(t[v].l, lo, mid) + value(t[v].r, mid, hi);
	}
};

int main() {
	rep(it,0,300) {
		int n = rand() % 7 + 1, m = rand() % 12, X = rand() % 3 ? 4 : 40;
		vector<vector<pii>> g(n);
		vector<vector<ll>> w(n, vector<ll>(n, -1));
		rep(i,0,m) {
			int a = rand() % n, b = rand() % n, x = rand() % X;
			g[a].push_back({b, x}); g[b].push_back({a, x});
			ll ww = 1LL << x;
			if (w[a][b] < 0 || ww < w[a][b]) w[a][b] = w[b][a] = ww;
		}
		int s = rand() % n;
		BigBits bb;
		auto r = dijkstra<int>(n, {{s, 0}},
			[&](int a, int b) { return bb.less(a, b); },
			[&](int u, int du, auto push) {
				for (auto [v, x] : g[u]) push(v, bb.add(du, x));
			});
		vector<vector<pair<int, ll>>> g2(n);
		rep(u,0,n) for (auto [v, x] : g[u]) g2[u].emplace_back(v, 1LL << x);
		auto d = dijkstra(g2, s);
		rep(i,0,n) {
			assert((d[i] != LLONG_MAX) == (bool)r.seen[i]);
			if (!r.seen[i]) continue;
			assert(bb.value(r.dist[i], 0, B) == d[i]);
			ll len = 0; int v = i;
			while (r.par[v] != -1) { assert(w[r.par[v]][v] >= 0); len += w[r.par[v]][v]; v = r.par[v]; }
			assert(v == s && len == d[i]);
		}
	}
	// carries across many bits: path of 1e5 edges of weight 2^0, plus one edge 2^16
	{
		int n = 100000;
		vector<vector<pii>> g(n);
		rep(i,0,n-1) g[i].push_back({i + 1, 0});
		g[0].push_back({n - 1, 16});
		BigBits bb;
		auto r = dijkstra<int>(n, {{0, 0}},
			[&](int a, int b) { return bb.less(a, b); },
			[&](int u, int du, auto push) {
				for (auto [v, x] : g[u]) push(v, bb.add(du, x));
			});
		assert(bb.value(r.dist[n-1], 0, B) == (1 << 16) && r.par[n-1] == 0);
		assert(bb.value(r.dist[n-2], 0, B) == n - 2);
	}
	cout << "Tests passed!" << endl;
}
