#include "../utilities/template.h"

#include "../../content/geometry/lineIntersection.h"
#include "../../content/geometry/LineProjRefl.h"
#include "../../content/geometry/lineDistance.h"
#include "../../content/geometry/SegmentDistance.h"

typedef long double ld;
typedef __int128 LL;
typedef Point<ll> PL;
typedef Point<double> PD;
typedef Point<ld> PE;

mt19937_64 rng(987654321);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }
PD toD(PL p) { return PD((double)p.x, (double)p.y); }
PE toE(PL p) { return PE((ld)p.x, (ld)p.y); }
bool inRange(PL p, ll C) { return abs(p.x) <= C && abs(p.y) <= C; }
LL abs128(LL x) { return x < 0 ? -x : x; }
int sgn128(LL x) { return (x > 0) - (x < 0); }

// r is num/den truncated toward zero, and nothing else
bool truncated(ll r, LL num, LL den) {
	if (den < 0) num = -num, den = -den;
	LL a = abs128(num), b = abs128((LL)r) * den;
	return b <= a && a < b + den && (r == 0 || sgn128(num) == sgn(r));
}

// exact intersection as (nx/det, ny/det), implicit line form
struct Exact { LL nx, ny, det; };
Exact interO(PL s1, PL e1, PL s2, PL e2) {
	LL a1 = e1.y - s1.y, b1 = s1.x - e1.x;
	LL c1 = a1 * s1.x + b1 * s1.y;
	LL a2 = e2.y - s2.y, b2 = s2.x - e2.x;
	LL c2 = a2 * s2.x + b2 * s2.y;
	return {c1 * b2 - c2 * b1, a1 * c2 - a2 * c1, a1 * b2 - a2 * b1};
}

ld segDistO(PL s, PL e, PL p) {
	LL vx = e.x - s.x, vy = e.y - s.y, wx = p.x - s.x, wy = p.y - s.y;
	LL ux = p.x - e.x, uy = p.y - e.y;
	LL vv = vx * vx + vy * vy, t = wx * vx + wy * vy;
	if (vv == 0 || t <= 0) return sqrtl((ld)(wx * wx + wy * wy));
	if (t >= vv) return sqrtl((ld)(ux * ux + uy * uy));
	return (ld)abs128(vx * wy - vy * wx) / sqrtl((ld)vv);
}

void lineInterDouble() {
	// reviewed reproducer: used to return (924262544, 1000000000)
	PD s1(999999997, 999999995), e1(999999996, 999999995);
	PD s2(975175172, 983029124), e2(975175175, 983029123);
	auto r = lineInter(s1, e1, s2, e2);
	assert(r.first == 1);
	assert(r.second.x == 924262559 && r.second.y == 999999995);
	auto q = lineInter(PE(s1.x, s1.y), PE(e1.x, e1.y),
		PE(s2.x, s2.y), PE(e2.x, e2.y));
	assert(q.first == 1 && q.second == PE(924262559, 999999995));

	// well conditioned lines far from the origin: the error must
	// be a few ulp of the coordinates, not of their cubes
	for (ll C : {0LL, 1000LL, 1000000LL, 1000000000LL}) {
		double worst = 0;
		rep(it,0,200000) {
			PL o(rnd(-C, C), rnd(-C, C));
			PL u(rnd(-10, 10), rnd(-10, 10)), w(rnd(-10, 10), rnd(-10, 10));
			ll cr = u.cross(w);
			if (cr == 0) continue;
			if ((ld)cr * cr * 100 < (ld)u.dist2() * w.dist2()) continue;
			PL a = o + PL(rnd(-1000, 1000), rnd(-1000, 1000));
			PL c = o + PL(rnd(-1000, 1000), rnd(-1000, 1000));
			PL b = a + u, d = c + w;
			Exact ex = interO(a, b, c, d);
			ld X = (ld)ex.nx / (ld)ex.det, Y = (ld)ex.ny / (ld)ex.det;
			auto res = lineInter(toD(a), toD(b), toD(c), toD(d));
			assert(res.first == 1);
			worst = max(worst, (double)max(abs(res.second.x - X),
				abs(res.second.y - Y)));
		}
		assert(worst <= 1e-9 + 4e-16 * (double)C);
	}
	// kinds
	assert(lineInter(PD(0,0), PD(1,1), PD(0,1), PD(1,2)).first == 0);
	assert(lineInter(PD(0,0), PD(1,1), PD(2,2), PD(5,5)).first == -1);
}

void lineInterInt() {
	const ll C = 1000000000;
	// integer answers at |coord| <= 1e9: must be exact
	// (products of three coordinates do not fit in long long)
	rep(it,0,400000) {
		ll m = it % 2 ? 30000 : 1000;
		PL x(rnd(-C, C), rnd(-C, C));
		PL u(rnd(-m, m), rnd(-m, m)), w(rnd(-m, m), rnd(-m, m));
		if (u.cross(w) == 0) continue;
		PL a = x + u * rnd(-m, m), b = x + u * rnd(-m, m);
		PL c = x + w * rnd(-m, m), d = x + w * rnd(-m, m);
		if (a == b || c == d) continue;
		if (!inRange(a, C) || !inRange(b, C)) continue;
		if (!inRange(c, C) || !inRange(d, C)) continue;
		auto res = lineInter(a, b, c, d);
		assert(res.first == 1 && res.second == x);
		assert(lineInter(c, d, b, a).second == x);
	}
	{ // corners of the range
		PL a(-C, -C), b(C, C), c(-C, C), d(C, -C);
		assert(lineInter(a, b, c, d) == make_pair(1, PL(0, 0)));
		// x = y = C/(2C-1) = 0.5000000003: truncated toward s1
		assert(lineInter(a, b, c, PL(C, -C + 2)) == make_pair(1, PL(0, 0)));
		assert(lineInter(b, a, c, PL(C, -C + 2)) == make_pair(1, PL(1, 1)));
		assert(lineInter(a, PL(-C, C), PL(C, C), c).second == PL(-C, C));
		assert(lineInter(a, b, PL(-C, -C + 1), PL(C, C + 1)).first == 0);
		assert(lineInter(a, b, PL(1, 1), PL(-7, -7)).first == -1);
	}
	// general answers: each coordinate is the exact one truncated
	// toward s1, whenever it fits
	for (ll R : {10LL, 1000000LL, C}) rep(it,0,300000) {
		PL a(rnd(-R, R), rnd(-R, R)), b(rnd(-R, R), rnd(-R, R));
		PL c(rnd(-R, R), rnd(-R, R)), d(rnd(-R, R), rnd(-R, R));
		if (it % 3 == 0) d = c + (b - a) + PL(rnd(-2, 2), rnd(-2, 2));
		if (a == b || c == d || !inRange(d, R)) continue;
		Exact ex = interO(a, b, c, d);
		auto res = lineInter(a, b, c, d);
		if (ex.det == 0) {
			assert(res.first == -((b - a).cross(c - a) == 0));
			continue;
		}
		assert(res.first == 1);
		LL lim = (LL)4000000000000000000LL * abs128(ex.det);
		if (abs128(ex.nx) > lim || abs128(ex.ny) > lim) continue;
		assert(truncated(res.second.x - a.x, ex.nx - a.x * ex.det, ex.det));
		assert(truncated(res.second.y - a.y, ex.ny - a.y * ex.det, ex.det));
	}
}

void lineProjTest() {
	const ll C = 1000000000;
	{
		PL a(-C, -C), b(C, C), p(C, -C);
		assert(lineProj(a, b, p) == PL(0, 0));
		assert(lineProj(a, b, p, true) == PL(-C, C));
		assert(lineProj(b, a, PL(-C, C), true) == p);
		assert(lineProj(a, PL(C, -C), PL(7, C)) == PL(7, -C));
		assert(lineProj(a, PL(C, -C), PL(7, C), true) == PL(7, -3*C));
	}
	rep(it,0,400000) { // integer answers at |coord| <= 1e9
		ll m = it % 2 ? 30000 : 100;
		PL a(rnd(-C, C), rnd(-C, C)), v(rnd(-m, m), rnd(-m, m));
		if (v == PL(0, 0)) continue;
		PL b = a + v * rnd(-m, m), f = a + v * rnd(-m, m);
		PL off = v.perp() * rnd(-m, m), p = f + off, r = f - off;
		if (a == b || !inRange(b, C) || !inRange(p, C)) continue;
		if (!inRange(r, C)) continue;
		assert(lineProj(a, b, p) == f);
		assert(lineProj(a, b, p, true) == r);
		assert(lineProj(b, a, r, true) == p);
	}
	// general answers: offset from p truncated toward zero
	for (ll R : {10LL, 1000000LL, C}) rep(it,0,300000) {
		PL a(rnd(-R, R), rnd(-R, R)), b(rnd(-R, R), rnd(-R, R));
		PL p(rnd(-R, R), rnd(-R, R));
		if (a == b) continue;
		LL vx = b.x - a.x, vy = b.y - a.y;
		LL cr = vx * (p.y - a.y) - vy * (p.x - a.x), d = vx * vx + vy * vy;
		rep(k,1,3) { // offset = -perp(v)*k*cr/d
			PL res = lineProj(a, b, p, k == 2);
			assert(truncated(res.x - p.x, vy * k * cr, d));
			assert(truncated(res.y - p.y, -vx * k * cr, d));
		}
	}
	// doubles far from the origin
	for (ll R : {1000LL, 1000000LL, C}) rep(it,0,200000) {
		PL o(rnd(-R, R), rnd(-R, R));
		PL a = o + PL(rnd(-1000, 1000), rnd(-1000, 1000));
		PL b = o + PL(rnd(-1000, 1000), rnd(-1000, 1000));
		PL p = o + PL(rnd(-1000, 1000), rnd(-1000, 1000));
		if (a == b) continue;
		PE v = toE(b - a), w = toE(p - a);
		PE f = toE(a) + v * (w.dot(v) / v.dist2());
		PD g = lineProj(toD(a), toD(b), toD(p));
		PD h = lineProj(toD(a), toD(b), toD(p), true);
		double tol = 1e-9 + 4e-16 * (double)R;
		assert(abs(g.x - f.x) <= tol && abs(g.y - f.y) <= tol);
		PE r = f * 2 - toE(p);
		assert(abs(h.x - r.x) <= tol && abs(h.y - r.y) <= tol);
	}
	PE e = lineProj(PE(0, 0), PE(3, 0), PE(1, 1) / 3);
	assert(e.x == 1 / 3.0L && e.y == 0);
}

void distTest() {
	static_assert(is_same<decltype(lineDist(PL(), PL(), PL())), double>::value, "");
	static_assert(is_same<decltype(lineDist(PD(), PD(), PD())), double>::value, "");
	static_assert(is_same<decltype(lineDist(PE(), PE(), PE())), ld>::value, "");
	static_assert(is_same<decltype(segDist(PL(), PL(), PL())), double>::value, "");
	static_assert(is_same<decltype(segDist(PD(), PD(), PD())), double>::value, "");
	static_assert(is_same<decltype(segDist(PE(), PE(), PE())), ld>::value, "");
	const ld eps = numeric_limits<ld>::epsilon();
	assert(abs(lineDist(PE(0, 0), PE(1, 1), PE(1, 0)) + sqrtl(0.5L)) <= eps);
	assert(abs(segDist(PE(0, 0), PE(1, 1), PE(1, 0)) - sqrtl(0.5L)) <= eps);
	assert(lineDist(PL(0, 0), PL(2, 0), PL(1, 1)) == 1);
	assert(lineDist(PL(0, 0), PL(2, 0), PL(1, -1)) == -1);

	{ // reviewed reproducer: relative error was 2.5e-4
		PL s(0, 0), e(933356, 835285), p(895611, 801506);
		ld o = segDistO(s, e, p);
		assert(abs(segDist(toD(s), toD(e), toD(p)) - o) <= 4e-16 * o);
		assert(abs(segDist(s, e, p) - o) <= 4e-16 * o);
	}
	{ // 1 unit off a segment of length 2.8e9
		const ll C = 1000000000;
		PL s(-C, -C), e(C - 1, C), p(0, 1);
		ld o = segDistO(s, e, p);
		assert(abs(o - 0.353553390328L) < 1e-11);
		assert(abs(segDist(s, e, p) - o) <= 4e-16 * o);
		assert(abs(lineDist(s, e, p) - o) <= 4e-16 * o);
		assert(abs(lineDist(e, s, p) + o) <= 4e-16 * o);
	}
	// no overflow from scaling by |e-s|^2 any more
	assert(segDist(PD(0, 0), PD(1e100, 0), PD(5e99, 1)) == 1);
	assert(segDist(PD(0, 0), PD(1e-100, 0), PD(5e-101, 1e-100)) == 1e-100);
	// degenerate and end cases
	assert(segDist(PL(5, 5), PL(5, 5), PL(8, 9)) == 5);
	assert(segDist(PL(0, 0), PL(2, 0), PL(5, 4)) == 5);
	assert(segDist(PL(0, 0), PL(2, 0), PL(-3, -4)) == 5);
	assert(segDist(PL(0, 0), PL(2, 0), PL(1, 0)) == 0);
	assert(segDist(PL(0, 0), PL(2, 0), PL(2, 0)) == 0);

	// long long input: exact up to the final division
	for (ll C : {5LL, 1000LL, 1000000LL, 1000000000LL}) rep(it,0,400000) {
		PL s(rnd(-C, C), rnd(-C, C)), e(rnd(-C, C), rnd(-C, C));
		PL p(rnd(-C, C), rnd(-C, C));
		if (it % 2 && !(s == e)) { // lattice point next to the segment
			PL v = e - s;
			ll t = rnd(0, 1000);
			p = s + PL((ll)((LL)v.x * t / 1000), (ll)((LL)v.y * t / 1000))
				+ PL(rnd(-2, 2), rnd(-2, 2));
			if (!inRange(p, C)) continue;
		}
		ld o = segDistO(s, e, p);
		double got = segDist(s, e, p);
		assert(abs(got - o) <= 4e-16 * o);
		assert(segDist(e, s, p) == got);
		if (!(s == e)) {
			double l = lineDist(s, e, p);
			assert(abs(l) <= got * (1 + 1e-15));
			assert(sgn(l) == sgn((e - s).cross(p - s)));
		}
		// doubles: exact while the cross product is (C <= 1e6),
		// absolute error ~ eps * C beyond that
		double gd = segDist(toD(s), toD(e), toD(p));
		if (C <= 1000000) assert(abs(gd - o) <= 4e-16 * o);
		else assert(abs(gd - o) <= 1e-15 * (double)C);
	}
}

// degenerate lines are points: -1 iff the point is on the other line
void degenerateLines() {
	typedef Point<ll> P;
	const int G = 3;
	rep(a,-G,G+1) rep(b,-G,G+1) rep(c,-G,G+1) rep(d,-G,G+1)
	rep(e,-G,G+1) rep(f,-G,G+1) {
		P p(a, b), s(c, d), t(e, f);
		int want = s == t ? -(p == s) : -(s.cross(t, p) == 0);
		assert(lineInter(p, p, s, t).first == want);
		assert(lineInter(s, t, p, p).first == want);
		Point<double> dp(a, b), ds(c, d), dt(e, f);
		assert(lineInter(dp, dp, ds, dt).first == want);
		assert(lineInter(ds, dt, dp, dp).first == want);
	}
	assert(lineInter(P(5,5), P(5,5), P(0,0), P(1,0)).first == 0);
}

int main() {
	degenerateLines();
	lineInterDouble();
	lineInterInt();
	lineProjTest();
	distTest();
	cout << "Tests passed!" << endl;
}
