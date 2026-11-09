#include "../utilities/template.h"

#include "../../content/data-structures/PersistentSegTree.h"

static mt19937_64 rng(1234567);
ll rnd(ll lo, ll hi) { return (ll)(rng() % (unsigned long long)(hi - lo + 1)) + lo; }

int main() {
	// query over branching versions, large values, clipped/empty/reversed ranges
	rep(it,0,4000) {
		int n = (int)rnd(1, 12);
		PST t(n);
		vector<vector<ll>> ver(1, vector<ll>(n, 0));
		vi roots{0};
		int ups = (int)rnd(0, 30);
		rep(u,0,ups) {
			int v = (int)rnd(0, sz(roots) - 1);
			int p = (int)rnd(0, n - 1);
			ll x = rnd(-(ll)4e17, (ll)4e17);
			if (rnd(0, 3) == 0) x = rnd(-3, 3);
			int before = sz(t.t);
			int r = t.upd(roots[v], p, x);
			int lg = 0; while ((1 << lg) < n) lg++;
			assert(sz(t.t) - before <= lg + 1);
			assert(r == before);
			roots.push_back(r);
			ver.push_back(ver[v]);
			ver.back()[p] += x;
		}
		rep(v,0,sz(roots)) {
			rep(b,-2,n+3) rep(e,-2,n+3) {
				ll s = 0;
				rep(i,max(b,0),min(e,n)) s += ver[v][i];
				assert(t.query(roots[v], b, e) == s);
			}
		}
	}

	// kth / counting on an array, vs sorting
	rep(it,0,4000) {
		int len = (int)rnd(0, 12), m = (int)rnd(1, 15);
		vi a(len);
		for (int& x : a) x = (int)rnd(0, m - 1);
		PST t(m); vi r{0};
		for (int x : a) r.push_back(t.upd(r.back(), x, 1));
		assert(sz(r) == len + 1);
		rep(lo,0,len+1) rep(hi,lo,len+1) {
			vi s(a.begin() + lo, a.begin() + hi);
			sort(all(s));
			rep(k,0,sz(s)) assert(t.kth(r[lo], r[hi], k) == s[k]);
			rep(x,0,m+1) rep(y,x,m+1) {
				ll c = 0;
				for (int v : s) c += x <= v && v < y;
				assert(t.query(r[hi], x, y) - t.query(r[lo], x, y) == c);
			}
		}
	}

	// weighted kth with huge multiplicities on branching versions
	rep(it,0,3000) {
		int n = (int)rnd(1, 10);
		PST t(n);
		vector<vector<ll>> ver(1, vector<ll>(n, 0));
		vi roots{0};
		// chain of versions so every pair (a, b) with a before b has b - a >= 0
		int ups = (int)rnd(1, 15);
		rep(u,0,ups) {
			int p = (int)rnd(0, n - 1);
			ll x = rnd(0, 3) == 0 ? rnd(1, 3) : rnd(1, (ll)1e17);
			roots.push_back(t.upd(roots.back(), p, x));
			ver.push_back(ver.back());
			ver.back()[p] += x;
		}
		rep(a,0,sz(roots)) rep(b,a,sz(roots)) {
			vector<ll> d(n);
			ll tot = 0;
			rep(i,0,n) tot += d[i] = ver[b][i] - ver[a][i];
			if (!tot) continue;
			rep(q,0,20) {
				ll k = rnd(0, tot - 1);
				if (q == 0) k = 0;
				if (q == 1) k = tot - 1;
				int p = 0; ll acc = d[0];
				while (acc <= k) acc += d[++p];
				assert(t.kth(roots[a], roots[b], k) == p);
			}
		}
	}

	// huge coordinate range: n = 2^30 and nearby, vs map oracle
	for (int n : {1 << 30, (1 << 30) - 1, 1000000007, 2, 1}) {
		PST t(n);
		vi roots{0};
		vector<map<int, ll>> ver(1);
		rep(u,0,200) {
			int p = rnd(0, 4) == 0 ? (rnd(0, 1) ? n - 1 : 0) : (int)rnd(0, n - 1);
			ll x = rnd(1, (ll)1e15);
			int before = sz(t.t);
			roots.push_back(t.upd(roots.back(), p, x));
			int lg = 0; while ((1LL << lg) < n) lg++;
			assert(sz(t.t) - before <= lg + 1);
			ver.push_back(ver.back());
			ver.back()[p] += x;
		}
		rep(q,0,3000) {
			int a = (int)rnd(0, sz(roots) - 1), b = (int)rnd(0, sz(roots) - 1);
			if (a > b) swap(a, b);
			int x = (int)rnd(0, n), y = (int)rnd(0, n);
			if (q % 7 == 0) x = 0;
			if (q % 11 == 0) y = n;
			ll s = 0;
			for (auto& [p, v] : ver[b]) if (x <= p && p < y) s += v;
			assert(t.query(roots[b], x, y) == s);
			ll tot = 0;
			for (auto& [p, v] : ver[b]) tot += v - ver[a][p];
			if (tot) {
				ll k = rnd(0, tot - 1);
				int pos = -1; ll acc = 0;
				for (auto& [p, v] : ver[b]) {
					acc += v - ver[a][p];
					if (acc > k) { pos = p; break; }
				}
				assert(t.kth(roots[a], roots[b], k) == pos);
			}
		}
	}

	// empty input: only root 0, everything is 0
	{
		PST t(5);
		assert(t.query(0, 0, 5) == 0 && t.query(0, -3, 9) == 0 && sz(t.t) == 1);
	}
	cout << "Tests passed!" << endl;
}
