#include "../utilities/template.h"

#include "../../content/geometry/CircleIntersection.h"
#include "../../content/geometry/CircleLine.h"
#include "../../content/geometry/CircleTangents.h"
#include "../../content/geometry/CirclePolyIntersection.h"
#include "../../content/geometry/MinimumEnclosingCircle.h"

// Regression tests for precision and degeneracy defects in the circle
// headers. Inputs are integers, decisions are checked against exact
// __int128 predicates and constructions against long double references
// computed in local coordinates.

typedef long double ld;
typedef __int128 lll;
typedef Point<ll> PL;

#define CHECK(c) do { if (!(c)) { fprintf(stderr, \
	"%s:%d: check failed: %s\n", __FILE__, __LINE__, #c); abort(); } } while (0)

mt19937_64 rng(20261012);
ll rnd(ll lo, ll hi) { return lo + (ll)(rng() % (unsigned long long)(hi - lo + 1)); }
P toP(PL p) { return P((double)p.x, (double)p.y); }
bool inRange(PL p, ll S) { return abs(p.x) <= S && abs(p.y) <= S; }
ld distL(P p, PL q) { return hypotl((ld)p.x - q.x, (ld)p.y - q.y); }
ld PIl = acosl(-1);

// Random integer vector (x, y) with integer length len <= maxLen.
PL pythag(ll maxLen, ll& len) {
	for (;;) {
		ll m = rnd(1, 60), n = rnd(0, m - 1);
		if (__gcd(m, n) != 1 || (m - n) % 2 == 0) continue;
		ll x = m*m - n*n, y = 2*m*n, l = m*m + n*n;
		if (l > maxLen) continue;
		ll k = rnd(1, maxLen / l);
		if (rng() % 4 == 0) k = maxLen / l;
		if (rng() % 2) swap(x, y);
		if (rng() % 2) x = -x;
		if (rng() % 2) y = -y;
		len = l * k;
		return PL(x * k, y * k);
	}
}

// Two integer circles within delta of external/internal tangency.
// Half of them have an integer centre distance, so delta = 0 is an exact
// tangency.
void nearTangentCircles(ll S, PL& a, PL& b, ll& r1, ll& r2) {
	for (;;) {
		PL d; ll len;
		if (rng() % 2) d = pythag(S, len);
		else {
			d = PL(rnd(-S, S), rnd(-S, S));
			len = llroundl(sqrtl((ld)d.dist2()));
		}
		a = PL(rnd(-S, S), rnd(-S, S)), b = a + d;
		if (len == 0 || !inRange(b, S)) continue;
		ll delta = rnd(-2, 2);
		if (rng() % 2) { // external: r1 + r2 = len + delta
			if (len < 2) continue;
			r1 = rnd(1, len - 1), r2 = len - r1 + delta;
		} else { // internal: r2 - r1 = len + delta
			if (len + 2 > S) continue;
			r1 = rnd(0, S - len - 2), r2 = r1 + len + delta;
		}
		if (r2 < 0 || r2 > S || r1 > S) continue;
		if (rng() % 2) swap(a, b), swap(r1, r2);
		return;
	}
}

void testPinned() {
	{ // exactly externally tangent, was reported as disjoint
		ll dx = -49962290LL - 252436090, dy = 594378458LL - 276860159;
		ll r1 = 68087418, r2 = 370390233;
		CHECK(dx*dx + dy*dy == (r1+r2)*(r1+r2));
		pair<P, P> out;
		CHECK(circleInter(P(252436090,276860159),
			P(-49962290,594378458), (double)r1, (double)r2, &out));
		CHECK((out.first - out.second).dist() < 1e-6);
		// one unit apart: disjoint
		CHECK(!circleInter(P(252436090,276860159),
			P(-49962290,594378458), (double)r1, (double)(r2-1), &out));
	}
	{ // exactly tangent line, was reported as missing the circle
		PL c(-173043,-54819), a(-98350,-37667), b(34223,88593);
		ll cr = (b-a).cross(c-a), len2 = (b-a).dist2(), r = 39092;
		CHECK((lll)cr*cr == (lll)r*r*len2);
		CHECK(sz(circleLine(toP(c), (double)r, toP(a), toP(b))) == 1);
		CHECK(sz(circleLine(toP(c), (double)(r-1), toP(a), toP(b))) == 0);
		CHECK(sz(circleLine(toP(c), (double)(r+1), toP(a), toP(b))) == 2);
	}
	{ // exactly internally tangent, was reported as nested
		ll dx = 581758624LL - 284653585, dy = 148560705LL + 247579347;
		ll r1 = 854864762, r2 = 359689697;
		CHECK(dx*dx + dy*dy == (r1-r2)*(r1-r2));
		P c1(284653585,-247579347), c2(581758624,148560705);
		CHECK(sz(tangents(c1, (double)r1, c2, (double)r2)) == 1);
		CHECK(sz(tangents(c1, (double)(r1-1), c2, (double)r2)) == 2);
		CHECK(sz(tangents(c1, (double)(r1+1), c2, (double)r2)) == 0);
	}
	{ // cross product exactly 1, radius was infinite
		double n = 5e8;
		P A(0, 0), B(n, n-1), C(n+1, n);
		ld ex = sqrtl((2*(ld)n*n - 2*n + 1) * (2*(ld)n*n + 2*n + 1) * 2) / 2;
		CHECK(abs(ccRadius(A, B, C) / (double)ex - 1) < 1e-14);
		P o = ccCenter(A, B, C);
		CHECK(abs((o - A).dist() / (double)ex - 1) < 1e-14);
		CHECK(abs((o - B).dist() / (double)ex - 1) < 1e-14);
		CHECK(abs((o - C).dist() / (double)ex - 1) < 1e-14);
	}
	{ // small circle centred on the edge of a huge square: was pi r^2
		vector<P> sq = {P(-1e9,-1e9), P(1e9,-1e9), P(1e9,1e9), P(-1e9,1e9)};
		for (double r : {1., 5., 20., 1000.}) {
			double half = (double)PIl * r * r / 2;
			CHECK(abs(circlePoly(P(1e9, 12345), r, sq) / half - 1) < 1e-9);
			CHECK(abs(circlePoly(P(1e9, 1e9), r, sq) / half - .5) < 1e-9);
			CHECK(abs(circlePoly(P(12345, -54321), r, sq) / half - 2) < 1e-9);
			CHECK(abs(circlePoly(P(1e9 + 2*r, 0), r, sq)) < 1e-9 * half);
		}
	}
	{ // duplicated point far from the origin: was NaN for some orders
		vector<P> ps;
		for (auto [x, y] : vector<pii>{{1,-1},{0,0},{1,0},{0,1},{1,-1},{-1,0}})
			ps.push_back(P(x + 1e9, y - 1e9));
		sort(all(ps));
		do {
			auto [o, r] = mec(ps);
			CHECK(abs(r - 1.1785113019775793) < 1e-12);
			CHECK((o - P(1e9 + 1./6, -1e9 - 1./6)).dist() < 1e-6);
		} while (next_permutation(all(ps)));
	}
	{ // a point 2.5e-9 r outside the circle through the other four
		vector<P> ps = {P(1e4,0), P(-1e4,0), P(0,1e4), P(0,-1e4), P(1e4,1)};
		rep(it,0,200) {
			auto [o, r] = mec(ps);
			CHECK(abs(r - 10000.00002499875) < 1e-9);
			for (P p : ps) CHECK((p - o).dist() < r * (1 + 1e-11));
		}
	}
}

void testCircleInter(ll S, int iters) {
	rep(it,0,iters) {
		PL a, b; ll r1, r2;
		nearTangentCircles(S, a, b, r1, r2);
		lll d2 = (lll)(b-a).dist2(), su = (lll)(r1+r2)*(r1+r2),
		    di = (lll)(r1-r2)*(r1-r2);
		bool expect = di <= d2 && d2 <= su;
		pair<P, P> out, out2;
		bool got = circleInter(toP(a), toP(b), (double)r1, (double)r2, &out);
		CHECK(got == expect);
		CHECK(circleInter(toP(b), toP(a), (double)r2, (double)r1, &out2) == expect);
		if (!got) continue;
		ld tol = 8e-16L * (ld)S + 1e-12L;
		for (P p : {out.first, out.second, out2.first, out2.second}) {
			CHECK(fabsl(distL(p, a) - r1) <= tol);
			CHECK(fabsl(distL(p, b) - r2) <= tol);
		}
		// .first is to the left of a->b when looking from a to b
		P ab = toP(b - a);
		CHECK(ab.cross(out.first - toP(a)) >= ab.cross(out.second - toP(a)));
		if (d2 == su || d2 == di) { // exact tangency: one known point
			ld len = sqrtl((ld)d2), k = (d2 == su || r1 > r2 ? r1 : -r1) / len;
			ld x = a.x + (b.x - a.x) * k, y = a.y + (b.y - a.y) * k;
			for (P p : {out.first, out.second})
				CHECK(hypotl(p.x - x, p.y - y) <= tol);
		}
	}
}

void testTangents(ll S, int iters) {
	rep(it,0,iters) {
		PL a, b; ll r1, r2;
		nearTangentCircles(S, a, b, r1, r2);
		for (ll s2 : {r2, -r2}) {
			lll h2 = (lll)(b-a).dist2() - (lll)(r1-s2)*(r1-s2);
			int expect = h2 < 0 ? 0 : h2 == 0 ? 1 : 2;
			auto tans = tangents(toP(a), (double)r1, toP(b), (double)s2);
			CHECK(sz(tans) == expect);
			ld tol = 8e-16L * (ld)S + 1e-12L;
			for (auto [p, q] : tans) {
				CHECK(fabsl(distL(p, a) - r1) <= tol);
				CHECK(fabsl(distL(q, b) - r2) <= tol);
				if (expect == 1) { CHECK((p - q).dist() <= 2 * tol); continue; }
				// the tangent is perpendicular to both radii, and both
				// circles lie on the side given by the signs of r1, r2
				ld len = hypotl((ld)q.x - p.x, (ld)q.y - p.y);
				ld ux = (q.x - p.x) / len, uy = (q.y - p.y) / len;
				ld s1 = ux * ((ld)a.y - p.y) - uy * ((ld)a.x - p.x);
				ld t2 = ux * ((ld)b.y - q.y) - uy * ((ld)b.x - q.x);
				// rounding of p, q is amplified by (centre distance) / len
				ld slack = tol * (1 + 2 * sqrtl((ld)(b-a).dist2()) / len);
				CHECK(fabsl(fabsl(s1) - r1) <= slack);
				CHECK(fabsl(fabsl(t2) - r2) <= slack);
				if (r1 > 0 && s2 != 0 && len > 1e-3L * S)
					CHECK((s1 > 0) == ((t2 > 0) == (s2 > 0)));
			}
			if (expect == 2 && h2 > (lll)S)
				CHECK(!(tans[0].first == tans[1].first) || r1 == 0);
		}
	}
}

void testCircleLine(ll S, int iters) {
	rep(it,0,iters) {
		PL a, b, c; ll r;
		if (rng() % 2) { // Pythagorean direction: exact tangency exists
			ll len; PL d = pythag(S, len);
			ll g = __gcd(abs(d.x), abs(d.y)), l0 = len / g;
			PL u(d.x / g, d.y / g);
			ll m = rng() % 2 ? rnd(-10, 10) : rnd(-S / l0, S / l0);
			if (m == 0) continue;
			a = PL(rnd(-S, S), rnd(-S, S)), b = a + d;
			c = a + u * rnd(-g, 2 * g) + u.perp() * m;
			r = abs(m) * l0 + rnd(-2, 2);
		} else {
			a = PL(rnd(-S, S), rnd(-S, S)), b = PL(rnd(-S, S), rnd(-S, S));
			c = PL(rnd(-S, S), rnd(-S, S));
			if (a == b) continue;
			r = llroundl(fabsl((ld)(b-a).cross(c-a)) / sqrtl((ld)(b-a).dist2()))
				+ rnd(-2, 2);
		}
		if (r < 0 || r > S || !inRange(b, S) || !inRange(c, S)) continue;
		ll s = (b-a).cross(c-a), L = (b-a).dist2();
		lll N = (lll)r*r*L - (lll)s*s;
		int expect = N < 0 ? 0 : N == 0 ? 1 : 2;
		auto res = circleLine(toP(c), (double)r, toP(a), toP(b));
		CHECK(sz(res) == expect);
		CHECK(sz(circleLine(toP(c), (double)r, toP(b), toP(a))) == expect);
		ld tol = 8e-16L * (ld)S + 1e-12L, len = sqrtl((ld)L);
		for (P p : res) {
			CHECK(fabsl(distL(p, c) - r) <= tol);
			ld off = ((ld)(b.x-a.x) * ((ld)p.y-a.y) - (ld)(b.y-a.y) * ((ld)p.x-a.x)) / len;
			CHECK(fabsl(off) <= tol);
		}
		if (expect == 2) { // ordered along a -> b
			CHECK(toP(b - a).dot(res[1] - res[0]) >= 0);
		}
		if (expect == 1) { // the projection of c
			ld t = (ld)(b-a).dot(c-a) / L;
			CHECK(hypotl(a.x + (b.x-a.x)*t - res[0].x, a.y + (b.y-a.y)*t - res[0].y) <= tol);
		}
	}
	// zero length line
	CHECK(circleLine(P(0, 0), 5, P(3, 4), P(3, 4)).empty());
	CHECK(circleLine(P(1e9, 1e9), 0, P(1e9, 1e9), P(1e9, 1e9)).empty());
}

void testCircumcircle(ll S, int iters) {
	rep(it,0,iters) {
		PL p[3];
		p[0] = PL(rnd(-S, S), rnd(-S, S));
		int mode = it % 3;
		if (mode == 0) rep(i,1,3) p[i] = PL(rnd(-S, S), rnd(-S, S));
		else if (mode == 1) { // cross product +-1
			ll n = rnd(1, S / 2);
			p[1] = p[0] + PL(n, n - 1), p[2] = p[0] + PL(n + 1, n);
			if (rng() % 2) rep(i,1,3) p[i] = p[0] + (p[i] - p[0]).perp();
		} else { // near collinear
			PL d(rnd(-S, S) / 3, rnd(-S, S) / 3);
			p[1] = p[0] + d, p[2] = p[0] + d * 2 + PL(rnd(-3, 3), rnd(-3, 3));
		}
		if (!inRange(p[1], S) || !inRange(p[2], S)) continue;
		PL b = p[2] - p[0], c = p[1] - p[0];
		ll cr = b.cross(c);
		if (cr == 0) continue;
		// offset of the centre from p[0], exact numerators
		lll nx = (lll)b.y * c.dist2() - (lll)c.y * b.dist2(),
		    ny = (lll)b.x * c.dist2() - (lll)c.x * b.dist2();
		ld ox = -(ld)nx / (2 * (ld)cr), oy = (ld)ny / (2 * (ld)cr);
		ld R = hypotl(ox, oy);
		int perm[3] = {0, 1, 2};
		do {
			P A = toP(p[perm[0]]), B = toP(p[perm[1]]), C = toP(p[perm[2]]);
			P o = ccCenter(A, B, C);
			double r = ccRadius(A, B, C);
			CHECK(fabsl(r / R - 1) < 1e-15L);
			CHECK(hypotl(o.x - (p[0].x + ox), o.y - (p[0].y + oy))
				<= 4e-16L * (R + S));
		} while (next_permutation(perm, perm + 3));
	}
}

// Area of the part of a disc of radius r at signed distance h inside a
// half-plane (h > 0: centre inside).
ld capArea(ld r, ld h) {
	if (h >= r) return PIl * r * r;
	if (h <= -r) return 0;
	return PIl * r * r - (r * r * acosl(h / r) - h * sqrtl(r * r - h * h));
}

void testCirclePoly(int iters) {
	const ll S = 1000000000;
	// Rotated and translated huge rectangles; the circle is small, so it
	// meets one edge (strip) or one corner (quadrant) only.
	rep(it,0,iters) {
		ll len; PL u = pythag(rng() % 2 ? 1 : 1000, len); // unit: u / len
		if (len == 1 && rng() % 2) u = PL(1, 0);
		ll W = rnd(S / 8, S / 4) / len, H = rnd(S / 8, S / 4) / len;
		PL t(rnd(-S / 4, S / 4), rnd(-S / 4, S / 4));
		auto tr = [&](ll x, ll y) { return t + u * x + u.perp() * y; };
		vector<P> poly = {toP(tr(-W, -H)), toP(tr(W, -H)), toP(tr(W, H)), toP(tr(-W, H))};
		rotate(poly.begin(), poly.begin() + rng() % 4, poly.end());
		ll rmax = rng() % 2 ? 10 : 100000;
		ll r = rnd(1, rmax), k = rnd(-r / len - 1, r / len + 1);
		int mode = (int)(rng() % 4);
		ld expect; PL c;
		if (mode == 0) { // near the edge y = -H, far from corners
			c = tr(rnd(-W / 2, W / 2), -H + k);
			expect = capArea((ld)r, (ld)k * len);
		} else if (mode == 1) { // near the edge x = W
			c = tr(W - k, rnd(-H / 2, H / 2));
			expect = capArea((ld)r, (ld)k * len);
		} else if (mode == 2) { // exactly at a corner
			c = tr(rng() % 2 ? W : -W, rng() % 2 ? H : -H);
			expect = PIl * r * r / 4;
		} else { // deep inside or far outside
			bool in = rng() % 2;
			c = tr(in ? rnd(-W / 2, W / 2) : W + 3 * r / len + 3, rnd(-H / 2, H / 2));
			expect = in ? PIl * r * r : 0;
		}
		double got = circlePoly(toP(c), (double)r, poly);
		CHECK(fabsl(got - expect) <= 1e-13L * PIl * r * r);
		reverse(all(poly)); // clockwise: negated
		CHECK(fabsl(circlePoly(toP(c), (double)r, poly) + expect)
			<= 1e-13L * PIl * r * r);
	}
	// Translating an integer instance must not change the answer, and
	// repeated vertices contribute nothing.
	rep(it,0,iters) {
		int n = (int)rnd(3, 8);
		vector<P> poly, moved, dup;
		P off((double)rnd(-S, S), (double)rnd(-S, S));
		rep(i,0,n) {
			poly.push_back(P((double)rnd(-20, 20), (double)rnd(-20, 20)));
			moved.push_back(poly.back() + off);
			dup.push_back(poly.back());
			if (rng() % 3 == 0) dup.push_back(poly.back());
		}
		P c((double)rnd(-20, 20), (double)rnd(-20, 20));
		double r = (double)rnd(0, 30), base = circlePoly(c, r, poly);
		CHECK(isfinite(base));
		CHECK(base == circlePoly(c + off, r, moved));
		CHECK(abs(base - circlePoly(c, r, dup)) < 1e-9);
	}
}

// Best of all circles through 1, 2 or 3 of the points that contain all
// points, in coordinates relative to p[0].
ld bruteMec(vector<PL> p) {
	int n = sz(p);
	ld best = 1e300L;
	PL p0 = p[0];
	for (PL& q : p) q = q - p0;
	auto ok = [&](ld x, ld y, ld r) {
		for (PL q : p)
			if (hypotl(q.x - x, q.y - y) > r * (1 + 1e-15L)) return false;
		return true;
	};
	rep(i,0,n) {
		if (ok((ld)p[i].x, (ld)p[i].y, 0)) best = 0;
		rep(j,0,i) {
			ld x = (p[i].x + p[j].x) / 2.0L, y = (p[i].y + p[j].y) / 2.0L;
			ld r = hypotl(p[i].x - x, p[i].y - y);
			if (r < best && ok(x, y, r)) best = r;
			rep(k,0,j) {
				PL b = p[j] - p[i], c = p[k] - p[i];
				ll cr = b.cross(c);
				if (cr == 0) continue;
				lll nx = (lll)b.y * c.dist2() - (lll)c.y * b.dist2(),
				    ny = (lll)b.x * c.dist2() - (lll)c.x * b.dist2();
				ld ox = -(ld)nx / (2 * (ld)cr), oy = (ld)ny / (2 * (ld)cr);
				r = hypotl(ox, oy);
				if (r < best && ok(p[i].x + ox, p[i].y + oy, r)) best = r;
			}
		}
	}
	return best;
}

void testMec(int iters) {
	CHECK(mec({}).second == 0);
	CHECK(mec({P(1e9, -1e9)}).second == 0);
	CHECK(mec({P(1e9, -1e9)}).first == P(1e9, -1e9));
	for (ll G : {1LL, 3LL, 10LL, 1000LL, 50000000LL}) rep(it,0,iters) {
		int n = (int)rnd(1, 8), mode = it % 5;
		vector<PL> p(n);
		rep(i,0,n) {
			p[i] = PL(rnd(-G, G), rnd(-G, G));
			if (mode == 1) p[i].y = p[i].x * 2; // collinear
			if (mode == 2 && i && rng() % 2) p[i] = p[rng() % i]; // duplicates
			if (mode == 3 && i > 1 && rng() % 2) // collinear triples
				p[i] = p[i-1] * 2 - p[i-2]; // |p[i]| <= 15 G
		}
		// any translation that keeps |coordinates| <= 1e9
		ll room = 1000000000 - 16 * G;
		PL off = mode == 4 ? PL(room, -room) : PL(rnd(-room, room), rnd(-room, room));
		if (it % 2 == 0) off = PL(0, 0);
		vector<P> ps;
		for (PL& q : p) q = q + off, ps.push_back(toP(q));
		ld ex = bruteMec(p);
		rep(rep_,0,3) { // mec shuffles internally: try several orders
			auto [o, r] = mec(ps);
			CHECK(isfinite(r) && isfinite(o.x) && isfinite(o.y));
			CHECK(fabsl(r - ex) <= 1e-11L * ex);
			// the centre is rounded to absolute doubles
			ld slack = 1e-11L * ex + 4e-16L * (fabsl((ld)off.x) + fabsl((ld)off.y) + 16 * G);
			for (PL q : p) CHECK(distL(o, q) <= r + slack);
			shuffle(all(ps), rng);
		}
	}
	// many identical and collinear points
	rep(it,0,200) {
		int n = (int)rnd(1, 300);
		PL base(rnd(-900000000, 900000000), rnd(-900000000, 900000000));
		PL d(rnd(-1000, 1000), rnd(-1000, 1000));
		vector<P> ps;
		ll lo = 0, hi = 0;
		rep(i,0,n) {
			ll k = it % 2 ? rnd(-5, 5) : rnd(-50000, 50000);
			if (i == 0) lo = hi = k;
			lo = min(lo, k), hi = max(hi, k);
			ps.push_back(toP(base + d * k));
		}
		auto [o, r] = mec(ps);
		ld ex = sqrtl((ld)d.dist2()) * (ld)(hi - lo) / 2;
		CHECK(isfinite(o.x) && isfinite(o.y));
		CHECK(fabsl(r - ex) <= 1e-11L * ex);
	}
}

int main() {
	testPinned();
	for (ll S : {1000LL, 1000000LL, 1000000000LL}) {
		testCircleInter(S, 100000);
		testTangents(S, 100000);
		testCircleLine(S, 100000);
		testCircumcircle(S, 50000);
	}
	testCirclePoly(100000);
	testMec(20000);
	// many tiny instances must stay cheap (no per-call setup cost)
	vector<P> tri = {P(0, 0), P(4, 0), P(0, 3)};
	double acc = 0;
	rep(it,0,200000) acc += mec(tri).second;
	CHECK(abs(acc - 200000 * 2.5) < 1e-6);
	cout<<"Tests passed!"<<endl;
}
