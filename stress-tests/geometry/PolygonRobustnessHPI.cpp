#include "../utilities/template.h"

#include "../../content/geometry/Point.h"
#include "../../content/geometry/PolygonArea.h"
#include "../../content/geometry/HalfPlaneIntersection.h"

// Exact brute force for small integer rays: the bounded intersection is
// the hull of all pairwise line intersections that satisfy every ray.
typedef long double ld;
typedef __int128 lll;
typedef Point<ll> IP;
typedef RayT<ll> IR;
mt19937_64 rng(20261012);
ll rnd(ll a, ll b) { return uniform_int_distribution<ll>(a, b)(rng); }

struct Fr { lll x, y, d; }; // (x/d, y/d), d > 0
int sg(lll x) { return (x > 0) - (x < 0); }
lll X(IP a, IP b) { return (lll)a.x * b.y - (lll)a.y * b.x; }
int orient(Fr a, Fr b, Fr c) {
	return sg((b.x*a.d - a.x*b.d) * (c.y*a.d - a.y*c.d)
		- (b.y*a.d - a.y*b.d) * (c.x*a.d - a.x*c.d));
}
bool frLess(Fr a, Fr b) {
	if (a.x*b.d != b.x*a.d) return a.x*b.d < b.x*a.d;
	return a.y*b.d < b.y*a.d;
}
// Expected vertices (CCW from the smallest), {} unless bounded with area.
vector<Fr> brute(const vector<IR>& rs) {
	int n = sz(rs);
	rep(i,0,n) for (int s = -1; s <= 1; s += 2) {
		bool allok = 1; // nonzero recession direction => unbounded or empty
		rep(j,0,n) if (X(rs[j].dp, rs[i].dp) * s < 0) allok = 0;
		if (allok) return {};
	}
	vector<Fr> v;
	rep(i,0,n) rep(j,i+1,n) {
		lll D = X(rs[i].dp, rs[j].dp);
		if (D == 0) continue;
		lll N = X(rs[j].p - rs[i].p, rs[j].dp);
		Fr f{rs[i].p.x * D + rs[i].dp.x * N, rs[i].p.y * D + rs[i].dp.y * N, D};
		if (f.d < 0) f.x = -f.x, f.y = -f.y, f.d = -f.d;
		bool in = 1;
		rep(k,0,n) if (rs[k].dp.x * (f.y - rs[k].p.y * f.d)
			- rs[k].dp.y * (f.x - rs[k].p.x * f.d) < 0) in = 0;
		if (in) v.push_back(f);
	}
	sort(all(v), frLess);
	v.erase(unique(all(v), [](Fr a, Fr b) {
		return !frLess(a, b) && !frLess(b, a); }), v.end());
	if (sz(v) < 3) return {};
	vector<Fr> h(sz(v)+1);
	int s = 0, t = 0;
	for (int it = 2; it--; s = --t, reverse(all(v)))
		for (Fr p : v) {
			while (t >= s + 2 && orient(h[t-2], h[t-1], p) <= 0) t--;
			h[t++] = p;
		}
	h.resize(t - (t == 2 && !frLess(h[0], h[1]) && !frLess(h[1], h[0])));
	if (sz(h) < 3) return {};
	return h;
}

template<class F>
bool finite(const vector<Point<F>>& v) {
	for (auto p : v) if (!isfinite(p.x) || !isfinite(p.y)) return false;
	return true;
}
// Same cyclic vertex sequence, within tol.
template<class F>
bool same(const vector<Point<F>>& got, const vector<Point<ld>>& want, ld tol) {
	if (sz(got) != sz(want)) return false;
	int n = sz(got);
	if (!n) return true;
	rep(s,0,n) {
		bool ok = 1;
		rep(i,0,n) {
			auto g = got[(i + s) % n];
			if (fabsl((ld)g.x - want[i].x) > tol || fabsl((ld)g.y - want[i].y) > tol) ok = 0;
		}
		if (ok) return true;
	}
	return false;
}
vector<Point<ld>> toLd(const vector<Fr>& h) {
	vector<Point<ld>> r;
	for (Fr f : h) r.push_back(Point<ld>((ld)f.x / (ld)f.d, (ld)f.y / (ld)f.d));
	return r;
}
vector<Ray> toDouble(const vector<IR>& rs, double dirScale = 1) {
	vector<Ray> r;
	for (auto& x : rs) r.push_back({Point<double>((double)x.p.x, (double)x.p.y),
		Point<double>((double)x.dp.x * dirScale, (double)x.dp.y * dirScale)});
	return r;
}
void dump(const vector<IR>& rs) {
	for (auto& r : rs) cerr << "(" << r.p.x << "," << r.p.y << ";" << r.dp.x << "," << r.dp.y << ") ";
	cerr << endl;
}

IP randPt(int R) { return IP(rnd(-R, R), rnd(-R, R)); }
IP randDir(int R) { IP d; do d = randPt(R); while (d == IP()); return d; }
// Random small instance: edges of a convex polygon and/or arbitrary rays.
vector<IR> gen(int R) {
	vector<IR> rs;
	int kind = (int)rnd(0, 3);
	if (kind) { // rays along the edges of a random convex polygon
		int k = (int)rnd(3, 7);
		vector<IP> pts;
		rep(i,0,k) pts.push_back(randPt(R));
		sort(all(pts));
		vector<IP> h(2*k);
		int s = 0, t = 0;
		for (int it = 2; it--; s = --t, reverse(all(pts)))
			for (IP p : pts) {
				while (t >= s + 2 && h[t-2].cross(h[t-1], p) <= 0) t--;
				h[t++] = p;
			}
		h.resize(max(t, 0));
		int m = sz(h);
		rep(i,0,m) if (m >= 2) {
			IP a = h[i], b = h[(i+1)%m];
			if (a == b) continue;
			IP dir = b - a;
			ll g = __gcd(abs(dir.x), abs(dir.y));
			if (rnd(0, 1)) dir = dir / g * rnd(1, 3);
			rs.push_back({rnd(0, 1) ? a : a + (b - a) / g * rnd(-2, 2), dir});
		}
	}
	int extra = kind == 1 ? 0 : (int)rnd(kind == 0 ? 1 : 0, kind == 0 ? 7 : 3);
	rep(i,0,extra) {
		if (sz(rs) && rnd(0, 3) == 0) { // parallel / antiparallel / duplicate
			IR r = rs[rnd(0, sz(rs)-1)];
			int w = (int)rnd(0, 3);
			if (w == 1) r.p = r.p + randPt(2);
			if (w == 2) r.dp = r.dp * -1;
			if (w == 3) r.dp = r.dp * -1, r.p = r.p + randPt(2);
			rs.push_back(r);
		} else rs.push_back({randPt(R), randDir(rnd(0, 2) ? 2 : R)});
	}
	shuffle(all(rs), rng);
	return rs;
}

ll cntPoly = 0, cntEmpty = 0, cntSliver = 0;
void checkSmall(const vector<IR>& rs) {
	auto want = toLd(brute(rs));
	(sz(want) ? cntPoly : cntEmpty)++;
	auto gi = halfPlaneIsect(rs);
	auto gd = halfPlaneIsect(toDouble(rs));
	if (!finite(gi) || !finite(gd) || !same(gi, want, 1e-12L) || !same(gd, want, 1e-9L)) {
		cerr << "small case mismatch: want " << sz(want) << " got " << sz(gi) << " / " << sz(gd) << endl;
		dump(rs);
		abort();
	}
	if (sz(gi)) assert(polygonArea2(gi) > 0);
	// scaling any direction, or everything, must not change the answer
	static const double scales[] = {1e-5, 1e-30, 1e5, 1e30};
	double s = scales[rnd(0, 3)];
	auto sd = toDouble(rs);
	for (auto& r : sd) if (rnd(0, 1)) r.dp = r.dp * s;
	if (!same(halfPlaneIsect(sd), want, 1e-9L)) {
		cerr << "direction scaling mismatch: want " << sz(want) << " got " << sz(halfPlaneIsect(sd)) << " scale " << s << endl;
		for (auto& r : sd) cerr << r.p << ";" << r.dp << " ";
		cerr << endl;
		abort();
	}
	double cs = scales[rnd(0, 3)];
	sd = toDouble(rs, s);
	for (auto& r : sd) r.p = r.p * cs;
	auto ws = want;
	for (auto& p : ws) p = p * (ld)cs;
	assert(same(halfPlaneIsect(sd), ws, 1e-9L * cs));
	// exact mode on a large integer image: p -> p*S + O, dp -> dp*S2
	ll S = rnd(1, 90000000), S2 = rnd(1, 200000000);
	IP O(rnd(-100000000, 100000000), rnd(-100000000, 100000000));
	auto big = rs;
	for (auto& r : big) r.p = r.p * S + O, r.dp = r.dp * S2;
	ws = want;
	for (auto& p : ws) p = p * (ld)S + Point<ld>((ld)O.x, (ld)O.y);
	auto gb = halfPlaneIsect(big);
	if (!same(gb, ws, 1e-6L)) {
		cerr << "scaled exact mismatch S=" << S << " S2=" << S2 << endl;
		dump(rs);
		abort();
	}
}

// Lattice triangle of area 1/2 with coordinates up to ~1e9.
void sliver() {
	ll a, b;
	do a = rnd(400000000, 999990000), b = rnd(400000000, 999990000); while (__gcd(a, b) != 1);
	// find c,d with a*d - b*c = 1
	ll x0 = 1, y0 = 0, x1 = 0, y1 = 1, r0 = a, r1 = b;
	while (r1) {
		ll q = r0 / r1;
		tie(r0, r1) = make_pair(r1, r0 - q * r1);
		tie(x0, x1) = make_pair(x1, x0 - q * x1);
		tie(y0, y1) = make_pair(y1, y0 - q * y1);
	}
	// a*x0 + b*y0 = 1  =>  d = x0, c = -y0
	ll c = -y0, d = x0;
	assert((lll)a * d - (lll)b * c == 1);
	IP O(rnd(-1000, 1000), rnd(-1000, 1000));
	// c,d have the same sign; keep all coordinates within 1e9
	if (c >= 0 && rnd(0, 1)) O = O - IP(a, b);
	vector<IP> tri = {O, O + IP(a, b), O + IP(c, d)};
	for (IP p : tri) if (abs(p.x) > 1000000000 || abs(p.y) > 1000000000) return;
	cntSliver++;
	assert(polygonArea2(tri) == 1);
	assert((polygonArea2<ll, lll>(tri)) == 1);
	vector<IR> ccw, cw;
	rep(i,0,3) {
		ccw.push_back({tri[i], tri[(i+1)%3] - tri[i]});
		cw.push_back({tri[(i+1)%3], tri[i] - tri[(i+1)%3]});
	}
	vector<Point<ld>> want;
	for (IP p : tri) want.push_back(Point<ld>((ld)p.x, (ld)p.y));
	rep(it,0,3) {
		shuffle(all(ccw), rng); shuffle(all(cw), rng);
		auto g = halfPlaneIsect(ccw);
		assert(same(g, want, 1e-6L));
		assert(fabsl(polygonArea2(g) - 1) < 1e-6L);
		assert(halfPlaneIsect(cw).empty());
	}
}

double area(const vector<Point<double>>& v) { return polygonArea2(v) / 2; }

int main() {
	typedef Point<double> P;
	// --- review reproducers ---
	// A3: unit square, direction vectors of any length
	for (double len : {1.0, 1e-4, 1e-5, 1e-9, 1e-100, 1e100}) {
		vector<Ray> sq = {{P(0,0),P(len,0)}, {P(1,0),P(0,len)}, {P(1,1),P(-len,0)}, {P(0,1),P(0,-len)}};
		auto r = halfPlaneIsect(sq);
		assert(sz(r) == 4 && abs(area(r) - 1) < 1e-12);
	}
	// B: contradictory parallel rays, no bounds
	assert(halfPlaneIsect(vector<Ray>{{P(0,1),P(1,0)}, {P(0,0),P(-1,0)}, {P(0,-1),P(-1,0)}}).empty());
	assert(halfPlaneIsect(vector<Ray>{{P(-1,1),P(1,0)}, {P(-2,-1),P(-1,0)}, {P(1,-1),P(2,0)}}).empty());
	assert(halfPlaneIsect(vector<IR>{{IP(0,1),IP(1,0)}, {IP(0,0),IP(-1,0)}, {IP(0,-1),IP(-1,0)}}).empty());
	// C: zero-area with add_bounds used to give NaN vertices
	{
		vector<Ray> c = {{P(4,2),P(-8,9)}, {P(0,8),P(0,9)}, {P(6,8),P(6,8)}, {P(3,10),P(-4,-7)}, {P(8,10),P(7,5)}};
		auto r = halfPlaneIsect(c, true);
		assert(r.empty());
	}
	// empty input, single rays, strips, unbounded wedges
	assert(halfPlaneIsect(vector<Ray>{}).empty());
	assert(halfPlaneIsect(vector<IR>{}).empty());
	assert(halfPlaneIsect(vector<Ray>{{P(0,0),P(1,0)}}).empty());
	assert(halfPlaneIsect(vector<Ray>{{P(0,0),P(1,0)}, {P(0,1),P(-1,0)}}).empty());
	assert(halfPlaneIsect(vector<Ray>{{P(0,0),P(1,0)}, {P(0,1),P(-1,0)}, {P(0,0),P(0,-1)}}).empty());
	assert(halfPlaneIsect(vector<Ray>{{P(0,0),P(1,0)}, {P(0,0),P(0,-1)}}).empty());
	{
		auto r = halfPlaneIsect(vector<Ray>{}, true);
		assert(sz(r) == 4 && abs(area(r) - 1e18) < 1e3);
		auto ri = halfPlaneIsect(vector<IR>{}, true);
		assert(sz(ri) == 4 && fabsl(polygonArea2(ri) - 2e18L) < 1);
		// strip 2 <= y <= 5 clipped by the bounds
		r = halfPlaneIsect(vector<Ray>{{P(0,2),P(3,0)}, {P(7,5),P(-2,0)}}, true);
		assert(sz(r) == 4 && abs(area(r) - 3e9) < 1e-3);
		// D: bounds only cover the non-negative quadrant
		r = halfPlaneIsect(vector<Ray>{{P(-1,-1),P(1,0)}, {P(1,-1),P(0,1)}, {P(1,1),P(-1,0)}, {P(-1,1),P(0,-1)}}, true);
		assert(sz(r) == 4 && abs(area(r) - 1) < 1e-12);
	}
	// a point and a segment are not polygons
	assert(halfPlaneIsect(vector<IR>{{IP(0,0),IP(1,0)}, {IP(0,0),IP(0,1)}, {IP(0,0),IP(-1,-1)}}).empty());
	assert(halfPlaneIsect(vector<IR>{{IP(0,0),IP(1,0)}, {IP(1,0),IP(0,1)}, {IP(1,0),IP(-1,0)}, {IP(0,0),IP(0,-1)}}).empty());
	// translated small polygon in double: local predicates keep it
	for (double M : {1e6, 1e9, 1e12}) {
		vector<Ray> sq = {{P(M,M),P(1,0)}, {P(M+1,M),P(0,1)}, {P(M+1,M+1),P(-1,0)}, {P(M,M+1),P(0,-1)}};
		auto r = halfPlaneIsect(sq);
		assert(sz(r) == 4 && abs(area(r) - 1) < 1e-9);
	}
	// extreme documented integer range (UBSan checks the __int128 products)
	rep(it,0,20000) {
		vector<IR> rs;
		int n = (int)rnd(1, 8);
		rep(i,0,n) {
			auto ext = [&](ll m) { return rnd(0, 2) ? (rnd(0, 1) ? m : -m) : rnd(-m, m); };
			IP dp;
			do dp = IP(ext(2000000000), ext(2000000000)); while (dp == IP());
			rs.push_back({IP(ext(1000000000), ext(1000000000)), dp});
		}
		assert(finite(halfPlaneIsect(rs)));
		assert(finite(halfPlaneIsect(rs, true)));
	}
	// --- brute force comparison ---
	rep(it,0,300000) checkSmall(gen((int)rnd(1, 5)));
	assert(cntPoly > 30000 && cntEmpty > 30000);
	rep(it,0,2000) sliver();
	assert(cntSliver > 1000);
	// --- many tiny calls and one large call stay fast ---
	{
		vector<IR> tri = {{IP(0,0),IP(1,0)}, {IP(1,0),IP(-1,1)}, {IP(0,1),IP(0,-1)}};
		ld tot = 0;
		rep(it,0,200000) { tri[0].p.y = it % 2 ? 0 : -1; tot += polygonArea2(halfPlaneIsect(tri)); }
		assert(tot > 0);
		int n = 200000;
		vector<IR> big;
		const ld PI = acosl(-1);
		rep(i,0,n) {
			ld a = 2 * PI * i / n;
			IP p((ll)llround(1e9 * cosl(a)), (ll)llround(1e9 * sinl(a)));
			ld b = a + PI / 2 + (ld)rnd(-1000, 1000) * 1e-7L;
			big.push_back({p, IP((ll)llround(1e9 * cosl(b)), (ll)llround(1e9 * sinl(b)))});
		}
		shuffle(all(big), rng);
		auto r = halfPlaneIsect(big);
		assert(sz(r) > 1000 && polygonArea2(r) > 6e18L && polygonArea2(r) < 6.3e18L);
	}
	cout<<"Tests passed!"<<endl;
}
