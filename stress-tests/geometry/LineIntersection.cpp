#include "../utilities/template.h"

#include "../../content/geometry/lineIntersection.h"
#include "../../content/geometry/lineDistance.h"

typedef long double ld;
typedef __int128 lll;
mt19937_64 rng(987654321);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }

// Oracle: write both lines as a*x + b*y = c and solve with Cramer's rule in
// 128-bit integers. Returns {type, numerators (X, Y), denominator D}.
struct Res { int type; lll X, Y, D; };
Res oracle(ll x1, ll y1, ll x2, ll y2, ll x3, ll y3, ll x4, ll y4) {
	lll a1 = y2 - y1, b1 = x1 - x2, c1 = a1 * x1 + b1 * y1;
	lll a2 = y4 - y3, b2 = x3 - x4, c2 = a2 * x3 + b2 * y3;
	lll D = a1 * b2 - a2 * b1;
	if (D == 0) {
		bool same = a1 * c2 == a2 * c1 && b1 * c2 == b2 * c1;
		return {same ? -1 : 0, 0, 0, 0};
	}
	return {1, c1 * b2 - c2 * b1, a1 * c2 - a2 * c1, D};
}

template<class T>
void test(ll x1, ll y1, ll x2, ll y2, ll x3, ll y3, ll x4, ll y4, ld tol) {
	if ((x1 == x2 && y1 == y2) || (x3 == x4 && y3 == y4)) return; // not lines
	typedef Point<T> P;
	Res r = oracle(x1, y1, x2, y2, x3, y3, x4, y4);
	auto got = lineInter(P((T)x1, (T)y1), P((T)x2, (T)y2), P((T)x3, (T)y3), P((T)x4, (T)y4));
	assert(got.first == r.type);
	if (r.type != 1) { assert(got.second == P(0, 0)); return; }
	if (T(0.5) != 0) { // floating point
		assert(fabsl((ld)got.second.x - (ld)r.X / (ld)r.D) <= tol);
		assert(fabsl((ld)got.second.y - (ld)r.Y / (ld)r.D) <= tol);
	} else if (r.X % r.D == 0 && r.Y % r.D == 0) { // integer intersection
		assert((lll)got.second.x == r.X / r.D);
		assert((lll)got.second.y == r.Y / r.D);
	}
}

int main() {
	rep(t,0,1000000) {
		const int GRID=10;
		Point<double>
			a(rand()%GRID, rand()%GRID),
			b(rand()%GRID, rand()%GRID),
			c(rand()%GRID, rand()%GRID),
			d(rand()%GRID, rand()%GRID);
		auto pa = lineInter(a,b,c,d);
		if (pa.first == 1) {
			assert(abs(lineDist(a, b, pa.second)) < 1e-8);
			assert(abs(lineDist(c, d, pa.second)) < 1e-8);
		}
	}
	// exhaustive tiny grid: parallel / coincident / unique all classified exactly
	const int G = 2;
	rep(x1,-G,G+1) rep(y1,-G,G+1) rep(x2,-G,G+1) rep(y2,-G,G+1)
	rep(x3,-G,G+1) rep(y3,-G,G+1) rep(x4,-G,G+1) rep(y4,-G,G+1) {
		test<double>(x1, y1, x2, y2, x3, y3, x4, y4, 1e-12);
		test<ll>(x1, y1, x2, y2, x3, y3, x4, y4, 0);
	}
	// random integer coordinates of growing magnitude
	for (ll r : {5LL, 100LL, 10000LL, 500000LL}) rep(it,0,500000) {
		ll x1 = rnd(-r, r), y1 = rnd(-r, r), x2 = rnd(-r, r), y2 = rnd(-r, r);
		ll x3 = rnd(-r, r), y3 = rnd(-r, r), x4 = rnd(-r, r), y4 = rnd(-r, r);
		int mode = it % 4;
		if (mode == 1) { // parallel or coincident
			ll dx = rnd(-r/4, r/4), dy = rnd(-r/4, r/4), k = rnd(-2, 2), m = rnd(-2, 2);
			x1 = rnd(-r/4, r/4), y1 = rnd(-r/4, r/4), x2 = x1 + dx, y2 = y1 + dy;
			x3 = rnd(0, 1) ? x1 + k * dx : rnd(-r/4, r/4);
			y3 = y1 + k * dy; x4 = x3 + m * dx, y4 = y3 + m * dy;
		}
		if (mode == 2) { // integer intersection point
			ll s = max(1LL, r / 8), X = rnd(-s, s), Y = rnd(-s, s);
			ll dx = rnd(-s, s), dy = rnd(-s, s), ex = rnd(-s, s), ey = rnd(-s, s);
			x1 = X + rnd(-3, 3) * dx, y1 = Y + (x1 - X) / (dx ? dx : 1) * dy;
			if (!dx) y1 = Y + rnd(-3, 3) * dy;
			x2 = x1 + dx, y2 = y1 + dy;
			x3 = X + rnd(-3, 3) * ex, y3 = Y + (x3 - X) / (ex ? ex : 1) * ey;
			if (!ex) y3 = Y + rnd(-3, 3) * ey;
			x4 = x3 + ex, y4 = y3 + ey;
		}
		// absolute error grows like 1/|sin(angle)|; scale tolerance by conditioning
		lll a1 = y2 - y1, b1 = x1 - x2, a2 = y4 - y3, b2 = x3 - x4, D = a1 * b2 - a2 * b1;
		ld cond = D ? sqrtl(ld(a1*a1 + b1*b1)) * sqrtl(ld(a2*a2 + b2*b2)) / fabsl((ld)D) : 1;
		test<double>(x1, y1, x2, y2, x3, y3, x4, y4, 1e-10 * (ld)r * cond);
		test<ll>(x1, y1, x2, y2, x3, y3, x4, y4, 0);
	}
	cout<<"Tests passed!"<<endl;
}
