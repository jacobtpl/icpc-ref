#include "../utilities/template.h"

#include "../../content/geometry/sideOf.h"

typedef long double ld;
mt19937_64 rng(4711);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }

// oracle: shoelace determinant of the triangle (s, e, p) in 128 bits
int orient(ll sx, ll sy, ll ex, ll ey, ll px, ll py) {
	__int128 d = (__int128)sx * ey - (__int128)ex * sy
		+ (__int128)ex * py - (__int128)px * ey
		+ (__int128)px * sy - (__int128)sx * py;
	return (d > 0) - (d < 0);
}
// signed distance from p to the line s-e
ld sdist(ld sx, ld sy, ld ex, ld ey, ld px, ld py) {
	ld a = ey - sy, b = sx - ex, c = -(a * sx + b * sy);
	return -(a * px + b * py + c) / hypotl(a, b);
}

int main() {
	typedef Point<ll> PL;
	typedef Point<double> PD;
	typedef Point<int> PI;
	// exhaustive small grid, all three coordinate types
	const int G = 3;
	rep(sx,-G,G+1) rep(sy,-G,G+1) rep(ex,-G,G+1) rep(ey,-G,G+1)
	rep(px,-G,G+1) rep(py,-G,G+1) {
		int o = orient(sx, sy, ex, ey, px, py);
		assert(sideOf(PL(sx,sy), PL(ex,ey), PL(px,py)) == o);
		assert(sideOf(PD(sx,sy), PD(ex,ey), PD(px,py)) == o);
		assert(sideOf(PI(sx,sy), PI(ex,ey), PI(px,py)) == o);
		// eps = 0 must agree with the exact version (also for s == e)
		assert(sideOf(PL(sx,sy), PL(ex,ey), PL(px,py), 0) == o);
		assert(sideOf(PD(sx,sy), PD(ex,ey), PD(px,py), 0) == o);
		if (sx == ex && sy == ey) continue;
		ld d = sdist(sx, sy, ex, ey, px, py);
		for (double eps : {0.1, 0.5, 0.75, 1.5, 3.3}) {
			if (fabsl(fabsl(d) - eps) < 1e-9) continue;
			int want = fabsl(d) <= eps ? 0 : d > 0 ? 1 : -1;
			assert(sideOf(PL(sx,sy), PL(ex,ey), PL(px,py), eps) == want);
			assert(sideOf(PD(sx,sy), PD(ex,ey), PD(px,py), eps) == want);
		}
	}
	// random long long points up to 1e9 (cross product < 2^63), with many
	// exactly and nearly collinear triples
	for (ll r : {5LL, 1000LL, 1000000LL, 1000000000LL}) rep(it,0,500000) {
		ll sx = rnd(-r, r), sy = rnd(-r, r), ex = rnd(-r, r), ey = rnd(-r, r);
		ll px = rnd(-r, r), py = rnd(-r, r);
		if (it % 3 == 0) {
			ll dx = rnd(-r/4, r/4), dy = rnd(-r/4, r/4), a = rnd(-2, 2), b = rnd(-2, 2);
			sx = rnd(-r/4, r/4), sy = rnd(-r/4, r/4);
			ex = sx + a * dx, ey = sy + a * dy;
			px = sx + b * dx + (it % 2 ? rnd(-1, 1) : 0), py = sy + b * dy;
		}
		int o = orient(sx, sy, ex, ey, px, py);
		assert(sideOf(PL(sx,sy), PL(ex,ey), PL(px,py)) == o);
		assert(sideOf(PL(sx,sy), PL(ex,ey), PL(px,py), 0) == o);
		assert(sideOf(PL(ex,ey), PL(sx,sy), PL(px,py)) == -o);
		if (r <= 1000000) { // exact in doubles
			assert(sideOf(PD((double)sx,(double)sy), PD((double)ex,(double)ey),
				PD((double)px,(double)py)) == o);
		}
		if (sx == ex && sy == ey) continue;
		ld d = sdist((ld)sx, (ld)sy, (ld)ex, (ld)ey, (ld)px, (ld)py);
		double eps = (double)rnd(0, 2 * r) / 7.0;
		if (fabsl(fabsl(d) - eps) < 1e-6 * (ld)r) continue;
		int want = fabsl(d) <= eps ? 0 : d > 0 ? 1 : -1;
		assert(sideOf(PL(sx,sy), PL(ex,ey), PL(px,py), eps) == want);
	}
	// random doubles with eps
	uniform_real_distribution<double> U(-100, 100), E(0, 50);
	rep(it,0,1000000) {
		double sx = U(rng), sy = U(rng), ex = U(rng), ey = U(rng), px = U(rng), py = U(rng);
		double eps = E(rng);
		ld d = sdist(sx, sy, ex, ey, px, py);
		if (fabsl(d) > 1e-9)
			assert(sideOf(PD(sx,sy), PD(ex,ey), PD(px,py)) == (d > 0 ? 1 : -1));
		if (fabsl(fabsl(d) - eps) < 1e-9) continue;
		int want = fabsl(d) <= eps ? 0 : d > 0 ? 1 : -1;
		assert(sideOf(PD(sx,sy), PD(ex,ey), PD(px,py), eps) == want);
	}
	// exactly at distance eps counts as on the line (3-4-5 triangles are exact)
	rep(k,1,1000) {
		assert(sideOf(PL(0,0), PL(3,4), PL(-4*k,3*k), 5*k) == 0);
		assert(sideOf(PL(0,0), PL(3,4), PL(4*k,-3*k), 5*k) == 0);
		assert(sideOf(PD(0,0), PD(3,4), PD(-4*k,3*k), 5*k) == 0);
		assert(sideOf(PD(0,0), PD(3,4), PD(-4*k-4,3*k+3), 5*k) == 1);
		assert(sideOf(PL(0,0), PL(3,4), PL(4*k+4,-3*k-3), 5*k) == -1);
		assert(sideOf(PL(1,1), PL(6,1), PL(k,3), 2) == 0);
		assert(sideOf(PL(1,1), PL(6,1), PL(k,4), 2) == 1);
	}
	cout<<"Tests passed!"<<endl;
}
