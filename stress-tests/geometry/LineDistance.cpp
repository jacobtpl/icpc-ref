#include "../utilities/template.h"

#include "../../content/geometry/lineDistance.h"

typedef long double ld;
mt19937_64 rng(2024);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }

// oracle: foot of the perpendicular in long double, sign from exact orientation
ld oracle(ll ax, ll ay, ll bx, ll by, ll px, ll py) {
	ld dx = ld(bx - ax), dy = ld(by - ay);
	ld t = (ld(px - ax) * dx + ld(py - ay) * dy) / (dx * dx + dy * dy);
	ld fx = ld(ax) + t * dx, fy = ld(ay) + t * dy;
	ld d = hypotl(ld(px) - fx, ld(py) - fy);
	__int128 o = (__int128)(bx - ax) * (py - ay) - (__int128)(by - ay) * (px - ax);
	return o > 0 ? d : o < 0 ? -d : 0;
}

template<class T>
void test(ll ax, ll ay, ll bx, ll by, ll px, ll py, ld scale) {
	typedef Point<T> P;
	P a((T)ax, (T)ay), b((T)bx, (T)by), p((T)px, (T)py);
	double got = lineDist(a, b, p);
	if (ax == bx && ay == by) { assert(std::isnan(got)); return; }
	ld want = oracle(ax, ay, bx, by, px, py);
	assert(fabsl(got - want) <= 1e-9 * scale);
	assert(sgn(got) == sgn(want));
	// antisymmetry
	assert(lineDist(b, a, p) == -got);
}

int main() {
	// exhaustive small grid
	const int G = 3;
	rep(ax,-G,G+1) rep(ay,-G,G+1) rep(bx,-G,G+1) rep(by,-G,G+1)
	rep(px,-G,G+1) rep(py,-G,G+1) {
		test<double>(ax, ay, bx, by, px, py, 1);
		test<ll>(ax, ay, bx, by, px, py, 1);
	}
	// random, increasing magnitude; ll products stay below 2^63 for |c| <= 1e9
	for (ll r : {10LL, 1000LL, 1000000LL, 1000000000LL}) rep(it,0,300000) {
		ll ax = rnd(-r, r), ay = rnd(-r, r), bx = rnd(-r, r), by = rnd(-r, r);
		ll px = rnd(-r, r), py = rnd(-r, r);
		if (it % 5 == 0) { // p (nearly) on the line
			ll k = rnd(-2, 2);
			px = ax + k * (bx - ax) + rnd(-1, 1), py = ay + k * (by - ay);
			if (max(abs(px), abs(py)) > r) continue;
		}
		test<ll>(ax, ay, bx, by, px, py, (ld)r);
		if (r <= 1000000) test<double>(ax, ay, bx, by, px, py, (ld)r);
	}
	// some fixed values
	typedef Point<double> P;
	assert(lineDist(P(0,0), P(2,0), P(1,3)) == 3);
	assert(lineDist(P(0,0), P(2,0), P(1,-3)) == -3);
	assert(lineDist(P(0,0), P(2,0), P(-100,0)) == 0);
	assert(abs(lineDist(P(0,0), P(1,1), P(0,1)) - sqrt(0.5)) < 1e-15);
	cout<<"Tests passed!"<<endl;
}
