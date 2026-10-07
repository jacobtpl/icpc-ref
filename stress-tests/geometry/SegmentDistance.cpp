#include "../utilities/template.h"

#include "../../content/geometry/SegmentDistance.h"

typedef long double ld;
mt19937_64 rng(3);
ll ri(ll a, ll b) { return uniform_int_distribution<ll>(a, b)(rng); }

// Oracle 1: clamp the projection parameter, in long double.
ld oracle(ld sx, ld sy, ld ex, ld ey, ld px, ld py) {
	ld dx = ex - sx, dy = ey - sy, l2 = dx*dx + dy*dy, t = 0;
	if (l2 > 0) t = min<ld>(1, max<ld>(0, ((px-sx)*dx + (py-sy)*dy) / l2));
	return hypotl(sx + t*dx - px, sy + t*dy - py);
}
// Oracle 2: ternary search over the segment (distance is convex in t).
ld ternary(P s, P e, P p) {
	ld lo = 0, hi = 1;
	auto f = [&](ld t) {
		return hypotl(s.x + t*((ld)e.x-s.x) - p.x, s.y + t*((ld)e.y-s.y) - p.y);
	};
	rep(it,0,200) {
		ld a = (2*lo + hi) / 3, b = (lo + 2*hi) / 3;
		if (f(a) < f(b)) hi = b; else lo = a;
	}
	return f(lo);
}

int main() {
	{ // usage example from the header
		Point<double> a, b(2,2), p(1,1);
		bool onSegment = segDist(a,b,p) < 1e-10;
		assert(onSegment);
	}
	// must be callable with temporaries / const points
	assert(abs(segDist(P(0,0), P(2,0), P(1,3)) - 3) < 1e-12);
	assert(abs(segDist(P(0,0), P(2,0), P(5,4)) - 5) < 1e-12);
	assert(abs(segDist(P(0,0), P(2,0), P(-3,-4)) - 5) < 1e-12);
	assert(abs(segDist(P(1,1), P(1,1), P(4,5)) - 5) < 1e-12);
	assert(segDist(P(1,1), P(1,1), P(1,1)) == 0);
	const P cs(0,0), ce(0,4), cp(3,2);
	assert(abs(segDist(cs, ce, cp) - 3) < 1e-12);

	// exhaustive on a small grid (lots of degenerate and collinear cases)
	const int G = 3;
	rep(a,-G,G+1) rep(b,-G,G+1) rep(c,-G,G+1) rep(d,-G,G+1) rep(x,-G,G+1) rep(y,-G,G+1) {
		P s(a,b), e(c,d), p(x,y);
		double r = segDist(s, e, p);
		assert(fabsl(r - oracle(a,b,c,d,x,y)) < 1e-12);
		if ((a + b + c + d + x + y) % 7 == 0)
			assert(fabsl(r - ternary(s, e, p)) < 1e-9);
	}
	// random integers of growing magnitude; the error scales with the coordinates
	for (ll lim : {10LL, 1000LL, 1000000LL}) rep(it,0,1000000) {
		P s((double)ri(-lim,lim), (double)ri(-lim,lim)), e((double)ri(-lim,lim), (double)ri(-lim,lim));
		P p((double)ri(-lim,lim), (double)ri(-lim,lim));
		if (it % 4 == 0) p = s + (e - s) * (double)ri(-2, 3); // on the line
		if (it % 9 == 0) e = s;
		double r = segDist(s, e, p);
		assert(fabsl(r - oracle(s.x,s.y,e.x,e.y,p.x,p.y)) < 1e-9 * (double)lim);
	}
	// random doubles
	uniform_real_distribution<double> U(-100, 100);
	rep(it,0,1000000) {
		P s(U(rng), U(rng)), e(U(rng), U(rng)), p(U(rng), U(rng));
		if (it % 5 == 0) e = s + (e - s) * 1e-6; // short segment
		double r = segDist(s, e, p);
		assert(fabsl(r - oracle(s.x,s.y,e.x,e.y,p.x,p.y)) < 1e-9);
		assert(r >= 0);
		assert(abs(r - segDist(e, s, p)) < 1e-9);
		if (it % 50 == 0) assert(fabsl(r - ternary(s, e, p)) < 1e-9);
	}
	cout<<"Tests passed!"<<endl;
}
