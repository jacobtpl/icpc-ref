#include "../utilities/template.h"

#include "../../content/geometry/circumcircle.h"

typedef long double ld;
mt19937_64 rng(5);
ll ri(ll a, ll b) { return uniform_int_distribution<ll>(a, b)(rng); }

// Oracle: solve the 2x2 linear system |X-A|^2 = |X-B|^2 = |X-C|^2 by
// Cramer's rule in long double, in absolute coordinates.
pair<ld, ld> oracle(P A, P B, P C) {
	ld a1 = 2*((ld)B.x-A.x), b1 = 2*((ld)B.y-A.y);
	ld c1 = (ld)B.x*B.x + (ld)B.y*B.y - (ld)A.x*A.x - (ld)A.y*A.y;
	ld a2 = 2*((ld)C.x-A.x), b2 = 2*((ld)C.y-A.y);
	ld c2 = (ld)C.x*C.x + (ld)C.y*C.y - (ld)A.x*A.x - (ld)A.y*A.y;
	ld det = a1*b2 - a2*b1;
	return {(c1*b2 - c2*b1) / det, (a1*c2 - a2*c1) / det};
}

void check(P A, P B, P C, double tol) {
	P c = ccCenter(A, B, C);
	double r = ccRadius(A, B, C);
	auto [ox, oy] = oracle(A, B, C);
	ld R = hypotl(ox - A.x, oy - A.y);
	ld err = max<ld>(1, R) * tol;
	assert(fabsl(c.x - ox) <= err && fabsl(c.y - oy) <= err);
	assert(fabsl(r - R) <= err);
	assert(fabsl((c - B).dist() - R) <= err && fabsl((c - C).dist() - R) <= err);
	// argument order must not matter
	P c2 = ccCenter(C, A, B);
	assert((c - c2).dist() <= err && abs(r - ccRadius(B, A, C)) <= err);
}

int main() {
	// known triangles; also checks temporaries/const arguments compile
	assert((ccCenter(P(0,0), P(4,0), P(0,3)) - P(2,1.5)).dist() < 1e-12);
	assert(abs(ccRadius(P(0,0), P(4,0), P(0,3)) - 2.5) < 1e-12);
	assert(abs(ccRadius(P(-1,0), P(1,0), P(0,sqrt(3))) - 2/sqrt(3)) < 1e-12);

	// all non-degenerate triangles on a small grid
	const int G = 3;
	rep(a,-G,G+1) rep(b,-G,G+1) rep(c,-G,G+1) rep(d,-G,G+1) rep(e,-G,G+1) rep(f,-G,G+1) {
		P A(a,b), B(c,d), C(e,f);
		if (A.cross(B, C) == 0) continue;
		check(A, B, C, 1e-12);
	}
	// random integers, including needle-like triangles and large coordinates
	for (ll lim : {10LL, 1000LL, 1000000LL}) rep(it,0,500000) {
		P A((double)ri(-lim,lim), (double)ri(-lim,lim)), B((double)ri(-lim,lim), (double)ri(-lim,lim));
		P C((double)ri(-lim,lim), (double)ri(-lim,lim));
		if (A.cross(B, C) == 0) continue;
		check(A, B, C, 1e-9);
	}
	// random doubles, well-conditioned (not too close to collinear)
	uniform_real_distribution<double> U(-100, 100);
	rep(it,0,500000) {
		P A(U(rng), U(rng)), B(U(rng), U(rng)), C(U(rng), U(rng));
		if (abs(A.cross(B, C)) < 1e-3) continue;
		check(A, B, C, 1e-8);
	}
	cout<<"Tests passed!"<<endl;
}
