#include "../utilities/template.h"

#include "../../content/geometry/PlanarFaces.h"
#include "../../content/geometry/SegmentIntersection.h"

typedef __int128 lll;

// Brute-force face tracer: for every half-edge, find the next one by an
// O(deg) atan2 scan, then walk cycles. Returns the set of walks (normalized).
set<vi> bruteFaces(const vector<P>& p, const vector<pii>& e) {
	int m = 2 * sz(e);
	auto from = [&](int h) { return h & 1 ? e[h/2].second : e[h/2].first; };
	auto to = [&](int h) { return h & 1 ? e[h/2].first : e[h/2].second; };
	vi nx(m);
	rep(h,0,m) {
		int v = to(h);
		double a = (p[from(h)] - p[v]).angle();
		double best = 1e18; int bh = -1;
		rep(k,0,m) if (from(k) == v && k != (h ^ 1)) {
			double b = (p[to(k)] - p[v]).angle();
			double d = fmod(a - b + 4 * M_PI, 2 * M_PI); // clockwise turn
			if (d < best) best = d, bh = k;
		}
		nx[h] = bh < 0 ? h ^ 1 : bh;
	}
	set<vi> res;
	vector<bool> seen(m);
	rep(i,0,m) if (!seen[i]) {
		vi w;
		for (int h = i; !seen[h]; h = nx[h]) seen[h] = 1, w.push_back(from(h));
		int b = int(min_element(all(w)) - w.begin());
		rotate(w.begin(), w.begin() + b, w.end());
		// walks may repeat a vertex; make the normal form canonical anyway
		vi best = w;
		rep(j,0,sz(w)) if (w[j] == w[0]) {
			vi c = w; rotate(c.begin(), c.begin() + j, c.end());
			best = min(best, c);
		}
		res.insert(best);
	}
	return res;
}

set<vi> normalize(const vector<vi>& f) {
	set<vi> res;
	for (vi w : f) {
		vi best; bool first = 1;
		rep(j,0,sz(w)) {
			vi c = w; rotate(c.begin(), c.begin() + j, c.end());
			if (first || c < best) best = c, first = 0;
		}
		res.insert(best);
	}
	return res;
}

// winding number of walk w about point q (q not on the walk)
int winding(const vector<P>& p, const vi& w, double qx, double qy) {
	int cnt = 0;
	rep(i,0,sz(w)) {
		P a = p[w[i]], b = p[w[(i + 1) % sz(w)]];
		double ax = (double)a.x - qx, ay = (double)a.y - qy;
		double bx = (double)b.x - qx, by = (double)b.y - qy;
		if ((ay > 0) != (by > 0)) {
			double cr = ax * by - ay * bx;
			if (ay <= 0 && by > 0 && cr > 0) cnt++;
			if (by <= 0 && ay > 0 && cr < 0) cnt--;
		}
	}
	return cnt;
}

lll shoelace(const vector<P>& p, const vi& w) {
	lll s = 0;
	rep(i,0,sz(w)) {
		P a = p[w[i]], b = p[w[(i + 1) % sz(w)]];
		s += (lll)a.x * b.y - (lll)a.y * b.x;
	}
	return s;
}

// checks everything that does not need small coordinates
void checkExact(const vector<P>& p, const vector<pii>& e, const PlanarFaces& F) {
	int n = sz(p), m = 2 * sz(e);
	assert(sz(F.fid) == m && sz(F.f) == sz(F.ar));
	// fid is consistent with the walks, every half-edge in exactly one walk
	map<pii, int> he;
	rep(i,0,sz(e)) he[e[i]] = 2 * i, he[pii(e[i].second, e[i].first)] = 2 * i + 1;
	vi cnt(m);
	rep(k,0,sz(F.f)) {
		const vi& w = F.f[k];
		assert(sz(w) >= 2);
		rep(j,0,sz(w)) {
			int h = he.at(pii(w[j], w[(j + 1) % sz(w)]));
			assert(F.fid[h] == k);
			cnt[h]++;
		}
		assert(shoelace(p, w) == (lll)F.ar[k]);
	}
	rep(h,0,m) assert(cnt[h] == 1);
	// components: Euler's formula, exactly one outer walk, ar == 0 iff tree
	vi comp(n, -1); int C = 0;
	vector<vi> adj(n);
	for (pii ed : e) adj[ed.first].push_back(ed.second), adj[ed.second].push_back(ed.first);
	rep(s,0,n) if (comp[s] < 0) {
		comp[s] = C; vi st = {s};
		while (!st.empty()) {
			int v = st.back(); st.pop_back();
			for (int u : adj[v]) if (comp[u] < 0) comp[u] = C, st.push_back(u);
		}
		C++;
	}
	vi V(C), E(C), W(C), outer(C), zero(C);
	rep(v,0,n) V[comp[v]]++;
	for (pii ed : e) E[comp[ed.first]]++;
	rep(k,0,sz(F.f)) {
		int c = comp[F.f[k][0]];
		W[c]++;
		if (F.ar[k] <= 0) outer[c]++;
		if (F.ar[k] == 0) zero[c]++;
	}
	rep(c,0,C) {
		if (E[c] == 0) { assert(W[c] == 0); continue; }
		assert(W[c] == E[c] - V[c] + 2);
		assert(outer[c] == 1);
		assert(zero[c] == (E[c] == V[c] - 1));
	}
}

// winding-number oracle; needs coordinates small enough for the 1e-3 offset
void checkGeometry(const vector<P>& p, const vector<pii>& e, const PlanarFaces& F) {
	int m = 2 * sz(e);
	vi comp(sz(p), -1); int C = 0;
	vector<vi> adj(sz(p));
	for (pii ed : e) adj[ed.first].push_back(ed.second), adj[ed.second].push_back(ed.first);
	rep(s,0,sz(p)) if (comp[s] < 0) {
		comp[s] = C; vi st = {s};
		while (!st.empty()) {
			int v = st.back(); st.pop_back();
			for (int u : adj[v]) if (comp[u] < 0) comp[u] = C, st.push_back(u);
		}
		C++;
	}
	rep(h,0,m) {
		P a = p[h & 1 ? e[h/2].second : e[h/2].first];
		P b = p[h & 1 ? e[h/2].first : e[h/2].second];
		P d = b - a;
		double len = d.dist();
		double qx = (double)(a.x + b.x) / 2 - 1e-3 * (double)d.y / len;
		double qy = (double)(a.y + b.y) / 2 + 1e-3 * (double)d.x / len;
		int k = F.fid[h];
		rep(j,0,sz(F.f)) {
			if (comp[F.f[j][0]] != comp[e[h/2].first]) continue;
			int wn = winding(p, F.f[j], qx, qy);
			if (F.ar[j] > 0) assert(wn == (j == k));
			else assert(wn == (F.ar[k] > 0 ? -1 : 0));
		}
	}
}

ll walks = 0, bounded = 0, multi = 0;
void checkAll(const vector<P>& p, const vector<pii>& e) {
	PlanarFaces F(p, e);
	walks += sz(F.f); for (ll a : F.ar) bounded += a > 0;
	multi += count_if(all(F.ar), [](ll a) { return a <= 0; }) > 1;
	checkExact(p, e, F);
	checkGeometry(p, e, F);
	assert(normalize(F.f) == bruteFaces(p, e));
	// same graph at extreme coordinates: identical walks, exact areas
	ll M = 1000000000, mx = 1;
	for (P q : p) mx = max({mx, abs(q.x), abs(q.y)});
	ll S = M / mx, O = M - mx * S;
	vector<ll> sx = {S, S, 1, S}, sy = {S, S, S, 1};
	vector<P> ox = {P(0,0), P(O,-O), P(-O,O), P(O,O)};
	rep(t,0,4) {
		vector<P> q(sz(p));
		rep(i,0,sz(p)) q[i] = P(p[i].x * sx[t] + ox[t].x, p[i].y * sy[t] + ox[t].y);
		PlanarFaces G(q, e);
		checkExact(q, e, G);
		assert(G.fid == F.fid && G.f == F.f);
		rep(k,0,sz(F.f)) assert((G.ar[k] > 0) == (F.ar[k] > 0) && (G.ar[k] == 0) == (F.ar[k] == 0));
	}
}

// random planar straight-line graph on the lattice [-G, G]^2
void randomCase(int n, int G, int tries) {
	n = min(n, (2*G+1) * (2*G+1));
	set<P> used; vector<P> p;
	while (sz(p) < n) {
		P q(rand() % (2*G+1) - G, rand() % (2*G+1) - G);
		if (used.insert(q).second) p.push_back(q);
	}
	vector<pii> e; set<pii> es;
	rep(it,0,tries) {
		int a = rand() % n, b = rand() % n;
		if (a == b || es.count(pii(a,b)) || es.count(pii(b,a))) continue;
		bool ok = 1;
		rep(v,0,n) if (v != a && v != b && onSegment(p[a], p[b], p[v])) ok = 0;
		for (pii ed : e) if (ok) {
			int c = ed.first, d = ed.second;
			if (c == a || c == b || d == a || d == b) continue;
			if (!segInter(p[a], p[b], p[c], p[d]).empty()) ok = 0;
		}
		if (!ok) continue;
		es.insert(pii(a,b));
		if (rand() & 1) swap(a, b);
		e.push_back(pii(a,b));
	}
	checkAll(p, e);
}

int main() {
	srand(2);
	// edge cases
	checkAll({}, {});
	checkAll({P(0,0)}, {});
	checkAll({P(0,0), P(1,1), P(-3,2)}, {});
	checkAll({P(0,0), P(1,0)}, {{0,1}});
	checkAll({P(0,0), P(1,0)}, {{1,0}});
	checkAll({P(0,0), P(0,1), P(0,2), P(0,3)}, {{0,1},{1,2},{2,3}});
	checkAll({P(0,0), P(1,1), P(2,2), P(3,3)}, {{3,2},{0,1},{1,2}});
	checkAll({P(0,0), P(2,0), P(0,2)}, {{0,1},{1,2},{2,0}});
	checkAll({P(0,0), P(2,0), P(0,2)}, {{1,0},{2,1},{0,2}});
	// triangle with a dangling edge inside and a component inside
	checkAll({P(0,0), P(9,0), P(0,9), P(1,1), P(2,2), P(3,1), P(4,1)},
	         {{0,1},{1,2},{2,0},{0,3},{3,4},{5,6}});
	// bow tie at a cut vertex
	checkAll({P(0,0), P(-2,-1), P(-2,1), P(2,-1), P(2,1)},
	         {{0,1},{1,2},{2,0},{0,3},{3,4},{4,0}});
	// all 8 compass directions + wheel
	{
		vector<P> p = {P(0,0)};
		rep(dx,-1,2) rep(dy,-1,2) if (dx || dy) p.push_back(P(dx, dy));
		vector<pii> e;
		rep(i,1,9) e.push_back({0, i});
		checkAll(p, e);
		vi ord = {1,2,3,5,8,7,6,4};
		rep(i,0,8) e.push_back({ord[i], ord[(i+1)%8]});
		checkAll(p, e);
	}
	{
		// nested squares and a diamond inside the inner one
		vector<P> p = {P(-9,-9),P(9,-9),P(9,9),P(-9,9), P(-5,-5),P(5,-5),P(5,5),P(-5,5),
		               P(0,-2),P(2,0),P(0,2),P(-2,0)};
		vector<pii> e = {{0,1},{1,2},{2,3},{3,0}, {4,5},{5,6},{6,7},{7,4}, {8,9},{9,10},{10,11},{11,8}};
		checkAll(p, e);
		e.push_back({0,4}); e.push_back({6,2});
		checkAll(p, e);
	}
	// extreme coordinates directly (not via checkAll's small-grid oracle)
	{
		ll M = 1000000000;
		vector<P> p = {P(-M,-M), P(M,-M), P(M,M), P(-M,M)};
		vector<pii> e = {{0,1},{1,2},{2,3},{3,0},{0,2}};
		PlanarFaces F(p, e);
		checkExact(p, e, F);
		assert(sz(F.f) == 3);
		ll tot = 0; int neg = 0;
		for (ll a : F.ar) { if (a < 0) neg++, assert(a == -8 * M * M); else tot += a; }
		assert(neg == 1 && tot == 8 * M * M);
		// square spiral path at +-1e9: a tree whose partial shoelace sums
		// overflow ll; the result must be exactly 0
		vector<P> q; vector<pii> f;
		ll r = M;
		q.push_back(P(-r,-r));
		rep(i,0,40) {
			int k = i % 4;
			ll s = r - (ll)(i / 4) * 50000000;
			ll s2 = r - (ll)((i + 1) / 4) * 50000000;
			if (k == 0) q.push_back(P(s, -s));
			if (k == 1) q.push_back(P(s, s));
			if (k == 2) q.push_back(P(-s, s));
			if (k == 3) q.push_back(P(-s2, -s2));
		}
		rep(i,0,sz(q)-1) f.push_back({i, i+1});
		PlanarFaces G(q, f);
		checkExact(q, f, G);
		assert(sz(G.f) == 1 && G.ar[0] == 0);
	}
	// random cases
	rep(it,0,20000) randomCase(rand() % 12 + 1, rand() % 5 + 1, rand() % 60);
	rep(it,0,3000) randomCase(rand() % 25 + 1, rand() % 3 + 1, 300);
	rep(it,0,2000) randomCase(rand() % 10 + 1, 30, 100);
	rep(it,0,200) randomCase(60, 6, 3000);
	cerr << walks << " walks, " << bounded << " bounded, " << multi << " multi-component cases" << endl;
	cout<<"Tests passed!"<<endl;
}
