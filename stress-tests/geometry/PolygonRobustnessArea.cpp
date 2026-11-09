#include "../utilities/template.h"

#include "../../content/geometry/Point.h"
#include "../../content/geometry/PolygonArea.h"
#include "../../content/geometry/PolygonCenter.h"
#include "../../content/geometry/PolygonCut.h"
#include "../../content/geometry/PolygonUnion.h"

// Translated / large-coordinate regressions for polygonArea2,
// polygonCenter, polyUnion and polygonCut, with exact integer oracles.
typedef long double ld;
typedef __int128 lll;
typedef Point<ll> IP;
typedef Point<double> DP;
mt19937_64 rng(20261014);
ll rnd(ll a, ll b) { return uniform_int_distribution<ll>(a, b)(rng); }

lll shoelace(const vector<IP>& v) { // origin based, exact
	lll a = 0;
	rep(i,0,sz(v)) {
		IP p = v[i], q = v[(i+1)%sz(v)];
		a += (lll)p.x * q.y - (lll)p.y * q.x;
	}
	return a;
}
vector<DP> toD(const vector<IP>& v, double ox = 0, double oy = 0) {
	vector<DP> r;
	for (IP p : v) r.push_back(DP((double)p.x + ox, (double)p.y + oy));
	return r;
}
// random simple polygon: star-shaped around the origin, CCW
vector<IP> star(int n, ll R) {
	for (;;) {
		vector<IP> v;
		rep(i,0,n) v.push_back(IP(rnd(-R, R), rnd(-R, R)));
		auto half = [](IP p) { return p.y > 0 || (p.y == 0 && p.x > 0); };
		bool bad = 0;
		for (IP p : v) if (p == IP()) bad = 1;
		sort(all(v), [&](IP a, IP b) {
			if (half(a) != half(b)) return half(a);
			lll c = (lll)a.x * b.y - (lll)a.y * b.x;
			return c != 0 ? c > 0 : a.dist2() < b.dist2(); });
		rep(i,0,n) { // consecutive rays must turn strictly left by < 180
			IP a = v[i], b = v[(i+1)%n];
			if ((lll)a.x * b.y - (lll)a.y * b.x <= 0) bad = 1;
		}
		if (!bad) return v;
	}
}
vector<IP> shift(vector<IP> v, IP o) { for (IP& p : v) p = p + o; return v; }

void testArea() {
	// review reproducer: unit square far from the origin
	for (double o : {0.0, 1e8, 1e9, 1e12, 4e15}) for (int sx : {-1, 1}) {
		double a = sx * o, b = o;
		vector<DP> sq = {DP(a,b), DP(a+1,b), DP(a+1,b+1), DP(a,b+1)};
		assert(polygonArea2(sq) == 2);
		reverse(all(sq));
		assert(polygonArea2(sq) == -2);
		rotate(sq.begin(), sq.begin() + 1, sq.end());
		assert(polygonArea2(sq) == -2);
	}
	assert(polygonArea2(vector<DP>{}) == 0);
	assert(polygonArea2(vector<IP>{}) == 0);
	assert(polygonArea2(vector<IP>{IP(3,4)}) == 0);
	assert(polygonArea2(vector<IP>{IP(3,4), IP(5,9)}) == 0);
	// exact integer path, small and at the documented limit
	ll M = (1LL << 30) - 1;
	{
		vector<IP> sq = {IP(-M,-M), IP(M,-M), IP(M,M), IP(-M,M)};
		assert(polygonArea2(sq) == 8 * M * M);
		assert((polygonArea2<ll, lll>(sq)) == (lll)8 * M * M);
		vector<IP> tri = {IP(M,M), IP(-M,M), IP(-M,-M)};
		assert(polygonArea2(tri) == 4 * M * M);
	}
	rep(it,0,200000) {
		int n = (int)rnd(3, 8);
		ll R = it % 2 ? 10 : 1000;
		vector<IP> v = star(n, R);
		lll want = shoelace(v);
		assert(want > 0);
		IP o = it % 4 < 2 ? IP() : IP(rnd(-M + R, M - R), rnd(-M + R, M - R));
		v = shift(v, o);
		rotate(v.begin(), v.begin() + rnd(0, n - 1), v.end());
		assert(polygonArea2(v) == (ll)want);
		assert((polygonArea2<ll, lll>(v)) == want);
		// small integers are exact in double at any representable offset
		double ox = (double)rnd(-1000000, 1000000) * 1e6, oy = (double)rnd(-1000000, 1000000) * 1e6;
		assert(polygonArea2(toD(v, ox, oy)) == (double)want);
		reverse(all(v));
		assert(polygonArea2(v) == -(ll)want);
	}
	// arbitrary (also self-intersecting) integer polygons, full range
	rep(it,0,200000) {
		int n = (int)rnd(0, 7);
		vector<IP> v;
		rep(i,0,n) v.push_back(IP(rnd(-M, M), rnd(-M, M)));
		assert((polygonArea2<ll, lll>(v)) == shoelace(v));
	}
}

// exact centroid as fractions X/(3A), Y/(3A) relative to the origin
void exactCenter(const vector<IP>& v, ld& cx, ld& cy) {
	lll A = 0, Xs = 0, Ys = 0;
	rep(i,0,sz(v)) {
		IP p = v[i], q = v[(i+1)%sz(v)];
		lll c = (lll)p.x * q.y - (lll)p.y * q.x;
		A += c; Xs += c * (p.x + q.x); Ys += c * (p.y + q.y);
	}
	cx = (ld)Xs / (ld)(3 * A); cy = (ld)Ys / (ld)(3 * A);
}
void testCenter() {
	// review reproducer: unit square far from the origin
	for (double o : {0.0, 1e8, 1e9, 1e12}) for (int sx : {-1, 1}) {
		double a = sx * o, b = o;
		vector<DP> sq = {DP(a,b), DP(a+1,b), DP(a+1,b+1), DP(a,b+1)};
		rep(r,0,4) {
			DP c = polygonCenter(sq);
			assert(c.x == a + 0.5 && c.y == b + 0.5);
			rotate(sq.begin(), sq.begin() + 1, sq.end());
		}
		reverse(all(sq));
		DP c = polygonCenter(sq);
		assert(c.x == a + 0.5 && c.y == b + 0.5);
	}
	{
		ll M = (1LL << 30) - 1;
		vector<IP> sq = {IP(-M,-M), IP(M,-M), IP(M,M), IP(-M,M)};
		auto c = polygonCenter(sq);
		assert(c.x == 0 && c.y == 0);
		vector<IP> tri = {IP(M,M), IP(M-3,M), IP(M,M-3)};
		c = polygonCenter(tri);
		assert(fabsl(c.x - (ld)(M-1)) < 1e-9L && fabsl(c.y - (ld)(M-1)) < 1e-9L);
		// zero area: documented NaN, no crash
		c = polygonCenter(vector<IP>{IP(0,0), IP(1,1), IP(2,2)});
		assert(isnan(c.x) && isnan(c.y));
		DP d = polygonCenter(vector<DP>{DP(0,0), DP(1,1), DP(2,2)});
		assert(isnan(d.x) && isnan(d.y));
	}
	ll M = (1LL << 30) - 1;
	rep(it,0,200000) {
		int n = (int)rnd(3, 8);
		ll R = it % 2 ? 10 : 1000;
		vector<IP> v = star(n, R);
		ld cx, cy;
		exactCenter(v, cx, cy); // small coordinates: no cancellation
		IP o = it % 4 < 2 ? IP() : IP(rnd(-M + R, M - R), rnd(-M + R, M - R));
		vector<IP> w = shift(v, o);
		rotate(w.begin(), w.begin() + rnd(0, n - 1), w.end());
		auto c = polygonCenter(w);
		assert(fabsl(c.x - (ld)o.x - cx) < 1e-8L && fabsl(c.y - (ld)o.y - cy) < 1e-8L);
		// double path: error relative to the polygon plus one rounding of the sum
		DP d = polygonCenter(toD(w));
		assert(fabsl(d.x - (ld)o.x - cx) < 1e-6L && fabsl(d.y - (ld)o.y - cy) < 1e-6L);
		reverse(all(w));
		c = polygonCenter(w);
		assert(fabsl(c.x - (ld)o.x - cx) < 1e-8L && fabsl(c.y - (ld)o.y - cy) < 1e-8L);
	}
}

// union area of small integer polygons by counting unit cells (rectilinear)
void testUnion() {
	// review reproducer: two overlapping unit squares far from the origin
	for (double o : {0.0, 1e8, 1e9, 1e12}) for (int sx : {-1, 1}) {
		double a = sx * o, b = o;
		vector<vector<DP>> ps = {
			{DP(a,b), DP(a+1,b), DP(a+1,b+1), DP(a,b+1)},
			{DP(a+0.5,b), DP(a+1.5,b), DP(a+1.5,b+1), DP(a+0.5,b+1)}};
		assert(polyUnion(ps) == 1.5);
		ps.push_back(ps[0]); // duplicate polygon
		assert(polyUnion(ps) == 1.5);
		ps.push_back({DP(a+5,b+5), DP(a+6,b+5), DP(a+5,b+6)});
		assert(polyUnion(ps) == 2);
	}
	{
		vector<vector<DP>> none;
		assert(polyUnion(none) == 0);
		vector<vector<DP>> emp = {{}, {DP(0,0), DP(1,0), DP(0,1)}, {}};
		assert(polyUnion(emp) == 0.5);
		vector<vector<IP>> iemp = {{}, {}};
		assert(polyUnion(iemp) == 0);
	}
	ll M = (1LL << 30) - 1;
	rep(it,0,30000) {
		// random axis-parallel rectangles and right triangles on a 12x12 grid
		const int G = 12;
		int n = (int)rnd(1, 5);
		vector<vector<IP>> ps;
		int twice[G][G] = {}; // covered area of each cell, times 2
		rep(k,0,n) {
			int x0 = (int)rnd(0, G-1), x1 = (int)rnd(x0+1, G), y0 = (int)rnd(0, G-1), y1 = (int)rnd(y0+1, G);
			if (rnd(0, 2)) {
				ps.push_back({IP(x0,y0), IP(x1,y0), IP(x1,y1), IP(x0,y1)});
				rep(x,x0,x1) rep(y,y0,y1) twice[x][y] = 2;
			} else { // staircase polygon: rectangle with a cut corner cell
				int cx = x1 - 1, cy = y1 - 1;
				if (x1 - x0 < 2 || y1 - y0 < 2) { k--; continue; }
				ps.push_back({IP(x0,y0), IP(x1,y0), IP(x1,cy), IP(cx,cy), IP(cx,y1), IP(x0,y1)});
				rep(x,x0,x1) rep(y,y0,y1) if (x != cx || y != cy) twice[x][y] = 2;
			}
		}
		int want2 = 0;
		rep(x,0,G) rep(y,0,G) want2 += twice[x][y];
		double want = want2 / 2.0;
		assert(abs(polyUnion(ps) - want) < 1e-9);
		// exact integer path anywhere in the documented range
		IP o(rnd(-M, M - G), rnd(-M, M - G));
		vector<vector<IP>> qs = ps;
		for (auto& p : qs) p = shift(p, o);
		assert(abs(polyUnion(qs) - want) < 1e-9);
		// doubles: integer offsets are represented exactly
		vector<vector<DP>> ds;
		for (auto& p : qs) ds.push_back(toD(p));
		assert(abs(polyUnion(ds) - want) < 1e-9);
		shuffle(all(qs), rng);
		for (auto& p : qs) rotate(p.begin(), p.begin() + rnd(0, sz(p) - 1), p.end());
		assert(abs(polyUnion(qs) - want) < 1e-9);
	}
	// general position: inclusion-exclusion of two convex polygons is not
	// available here, so check translation invariance on random triangles
	rep(it,0,30000) {
		int n = (int)rnd(1, 4);
		vector<vector<IP>> ps;
		rep(k,0,n) ps.push_back(star(3, 20));
		double base = polyUnion(ps);
		IP o(rnd(-M + 20, M - 20), rnd(-M + 20, M - 20));
		vector<vector<IP>> qs = ps;
		for (auto& p : qs) p = shift(p, o);
		assert(abs(polyUnion(qs) - base) < 1e-9);
		vector<vector<DP>> ds;
		for (auto& p : qs) ds.push_back(toD(p));
		assert(abs(polyUnion(ds) - base) < 1e-9);
		if (n == 1) assert(abs(base - (double)shoelace(ps[0]) / 2) < 1e-9);
	}
}

void testCut() {
	// thin translated polygons: compare each new vertex with the exact
	// rational intersection of the integer inputs
	ll cnt = 0;
	rep(it,0,200000) {
		ll O = it % 3 == 0 ? 0 : it % 3 == 1 ? 10000000 : 100000000;
		IP o(rnd(-O, O), rnd(-O, O));
		ll L = rnd(1, 1000);
		// thin triangle/quad of width ~1 and length L
		vector<IP> v = {o, o + IP(L, rnd(0, 1)), o + IP(L, 2), o + IP(rnd(0, 1), 1)};
		if (shoelace(v) <= 0) continue;
		IP s = o + IP(rnd(-5, L + 5), rnd(-5, 5)), e = o + IP(rnd(-5, L + 5), rnd(-5, 5));
		if (s == e) continue;
		auto side = [&](IP p) {
			return (lll)(e.x - s.x) * (p.y - s.y) - (lll)(e.y - s.y) * (p.x - s.x); };
		vector<pair<ld, ld>> want;
		rep(i,0,4) {
			IP cur = v[i], prev = v[(i+3)%4];
			lll a = side(cur), b = side(prev);
			if ((a < 0) != (b < 0)) {
				// cur + (prev - cur) * a / (a - b), relative to o
				ld t = (ld)a / (ld)(a - b);
				want.emplace_back((ld)(cur.x - o.x) + (ld)(prev.x - cur.x) * t,
					(ld)(cur.y - o.y) + (ld)(prev.y - cur.y) * t);
			}
			if (a < 0) want.emplace_back((ld)(cur.x - o.x), (ld)(cur.y - o.y));
		}
		auto got = polygonCut(toD(v), DP((double)s.x, (double)s.y), DP((double)e.x, (double)e.y));
		assert(sz(got) == sz(want));
		rep(i,0,sz(got)) {
			// one ulp of the offset (1.5e-8 at 1e8) plus local rounding
			assert(fabsl((ld)got[i].x - (ld)o.x - want[i].first) < 1e-7L);
			assert(fabsl((ld)got[i].y - (ld)o.y - want[i].second) < 1e-7L);
		}
		cnt += sz(want) > 0;
		// the cut part's area must agree too (local accumulation)
		vector<pair<ld, ld>> w = want;
		ld ar = 0;
		rep(i,1,sz(w)-1) ar += (w[i].first - w[0].first) * (w[i+1].second - w[0].second)
			- (w[i].second - w[0].second) * (w[i+1].first - w[0].first);
		assert(fabsl((ld)polygonArea2(got) - ar) < 1e-4L);
	}
	assert(cnt > 50000);
	assert(polygonCut({}, DP(0,0), DP(1,0)).empty());
}

int main() {
	testArea();
	testCenter();
	testUnion();
	testCut();
	// tiny repeat
	{
		vector<DP> sq = {DP(0,0), DP(1,0), DP(1,1), DP(0,1)};
		vector<vector<DP>> two = {sq, sq};
		double s = 0;
		rep(it,0,300000) {
			s += polygonArea2(sq) + polygonCenter(sq).x + polyUnion(two)
				+ sz(polygonCut(sq, DP(0.5,0), DP(0.5,1)));
		}
		assert(s == 300000 * (2 + 0.5 + 1 + 4.0));
	}
	// large inputs
	{
		int n = 1000000;
		vector<IP> v;
		rep(i,0,n) { ll x = i - n / 2; v.push_back(IP(x, x * x / 1000 - 500000000)); }
		assert((polygonArea2<ll, lll>(v)) == shoelace(v));
		auto c = polygonCenter(v);
		ld cx, cy; exactCenter(v, cx, cy);
		assert(fabsl(c.x - cx) < 1e-6L && fabsl(c.y - cy) < 1e-6L);
		vector<vector<IP>> ps;
		rep(i,0,300) ps.push_back(shift(star(4, 30), IP(rnd(-100, 100) + 1000000000, rnd(-100, 100))));
		double u = polyUnion(ps);
		assert(u > 0 && u < 260.0 * 260.0);
	}
	cout << "Tests passed!" << endl;
}
