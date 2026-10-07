#include "../utilities/template.h"

#include "../../content/geometry/PolygonArea.h"
#include "../../content/geometry/PolygonCenter.h"
#include "../../content/geometry/InsidePolygon.h"

namespace hull {
#include "../../content/geometry/ConvexHull.h"
}

// Pick's theorem on a small lattice polygon: 2A = 2I + B - 2, with I and B
// counted by brute force over the bounding box.
void testPick() {
	typedef Point<ll> P;
	mt19937 rng(4);
	auto ri = [&](int a, int b) { return uniform_int_distribution<int>(a, b)(rng); };
	int nondeg = 0;
	rep(it,0,100000) {
		int C = ri(1, 7), N = ri(1, 10);
		vector<P> ps;
		rep(i,0,N) ps.emplace_back(ri(-C, C), ri(-C, C));
		vector<P> h = hull::convexHull(ps);
		rotate(h.begin(), h.begin() + ri(0, sz(h) - 1), h.end());
		int n = sz(h);
		ll want = 0;
		if (n >= 3) {
			ll I = 0, B = 0;
			rep(x,-C,C+1) rep(y,-C,C+1) {
				bool in = true, on = false;
				rep(i,0,n) {
					ll c = h[i].cross(h[(i + 1) % n], P(x, y));
					if (c < 0) in = false;
					if (c == 0) on = true;
				}
				if (in && on) B++;
				else if (in) I++;
			}
			want = 2 * I + B - 2;
			nondeg++;
		}
		assert(polygonArea2(h) == want);
		reverse(all(h));
		assert(polygonArea2(h) == -want); // clockwise
		// translation invariance with large offsets
		P off(ri(-1000000000, 1000000000), ri(-1000000000, 1000000000));
		for (P& p : h) p = p + off;
		assert(polygonArea2(h) == -want);
	}
	assert(nondeg > 50000);
}

// non-convex (even self-intersecting) closed polylines: compare with the
// trapezoid formula in 128 bit, |x|,|y| <= 1e9
void testTrapezoid() {
	typedef Point<ll> P;
	typedef __int128 L;
	mt19937_64 rng(6);
	auto ri = [&](ll a, ll b) { return uniform_int_distribution<ll>(a, b)(rng); };
	for (ll M : {3LL, 1000LL, 1000000000LL}) rep(it,0,200000) {
		int n = (int)ri(1, 12);
		vector<P> v;
		rep(i,0,n) v.emplace_back(ri(-M, M), ri(-M, M));
		L want = 0;
		rep(i,0,n) {
			P a = v[i], b = v[(i + 1) % n];
			want += (L)(a.x - b.x) * (a.y + b.y);
		}
		assert((L)polygonArea2(v) == want);
	}
	// axis-aligned rectangle at the coordinate limit and a big staircase
	vector<P> sq = {P(-1000000000, -1000000000), P(1000000000, -1000000000),
		P(1000000000, 1000000000), P(-1000000000, 1000000000)};
	assert(polygonArea2(sq) == 8000000000000000000LL);
	const int K = 100000;
	vector<P> st;
	ll area = 0;
	rep(i,0,K) st.emplace_back(i, i), st.emplace_back(i + 1, i), area += i + 1;
	st.emplace_back(K, K), st.emplace_back(0, K);
	area = (ll)K * K - area + K; // square minus the part below the staircase
	assert(polygonArea2(st) == 2 * area);
	// doubles
	vector<Point<double>> tri = {Point<double>(0.5, 0.5), Point<double>(2.5, 0.5), Point<double>(0.5, 1.5)};
	assert(polygonArea2(tri) == 2.0);
}

int main() {
	testPick();
	testTrapezoid();
	srand(0);
	typedef Point<double> P;
	vector<P> ps = {P{0,0}, P{6,4}, P{0,9}};
	int count = 0;
	P su{0,0};
	rep(it,0,100000) {
		double x = rand() / (RAND_MAX + 1.0);
		double y = rand() / (RAND_MAX + 1.0);
		x *= 10;
		y *= 10;
		if (!inPolygon(ps, P{x,y}, true)) continue;
		count++;
		su = su + P{x,y};
	}
	su = su / count;
	double approxArea = (double)count / 100000 * 100;
	assert(abs(polygonArea2(ps)/2.0 - approxArea) < 1);
	auto p = polygonCenter(ps);
	assert(abs(p.x - su.x) < 1e-1 && abs(p.y - su.y) < 1e-1);
	cout<<"Tests passed!"<<endl;
}
