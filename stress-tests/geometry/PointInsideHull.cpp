#include "../utilities/template.h"

#include "../../content/geometry/ConvexHull.h"
namespace ignore {
	#include "../../content/geometry/SegmentDistance.h"
}
#include "../../content/geometry/PointInsideHull.h"
#include "../../content/geometry/InsidePolygon.h"

// -1 outside, 0 on the boundary, 1 strictly inside; O(n), for a strictly
// convex CCW polygon (or a point / segment when n < 3)
int slow(const vector<P>& l, P p) {
	int n = sz(l);
	typedef __int128 L;
	if (n < 3) {
		P s = l[0], e = l.back();
		L c = (L)(e.x - s.x) * (p.y - s.y) - (L)(e.y - s.y) * (p.x - s.x);
		return c == 0 && min(s.x, e.x) <= p.x && p.x <= max(s.x, e.x)
			&& min(s.y, e.y) <= p.y && p.y <= max(s.y, e.y) ? 0 : -1;
	}
	int res = 1;
	rep(i,0,n) {
		P a = l[i], b = l[(i + 1) % n];
		L c = (L)(b.x - a.x) * (p.y - a.y) - (L)(b.y - a.y) * (p.x - a.x);
		if (c < 0) return -1;
		if (c == 0) res = 0;
	}
	return res;
}

void checkHull(const vector<P>& h, P p) {
	int s = slow(h, p);
	assert(inHull(h, p, true) == (s > 0));
	assert(inHull(h, p, false) == (s >= 0));
}

void testAgainstBrute() {
	mt19937_64 rng(11);
	auto ri = [&](ll a, ll b) { return uniform_int_distribution<ll>(a, b)(rng); };
	// tiny hulls (incl. single points and segments), every start vertex, every grid point
	rep(it,0,30000) {
		int C = (int)ri(1, 5), N = (int)ri(1, 10);
		vector<P> ps;
		rep(i,0,N) ps.emplace_back(ri(-C, C), ri(-C, C));
		vector<P> h = convexHull(ps);
		rotate(h.begin(), h.begin() + ri(0, sz(h) - 1), h.end());
		rep(x,-C-1,C+2) rep(y,-C-1,C+2) checkHull(h, P(x, y));
	}
	// large coordinates (|x|,|y| <= 1e9); query vertices, edge points and random points
	const ll M = 1000000000;
	rep(it,0,100000) {
		int N = (int)ri(1, 30);
		vector<P> ps;
		rep(i,0,N) ps.emplace_back(ri(-M, M), ri(-M, M));
		if (it % 4 == 0) for (auto &p : ps) p = P(p.x / 2 * 2, p.y / 2 * 2); // even: midpoints are lattice points
		vector<P> h = convexHull(ps);
		rotate(h.begin(), h.begin() + ri(0, sz(h) - 1), h.end());
		int n = sz(h);
		rep(q,0,10) {
			int i = (int)ri(0, n - 1), j = (i + 1) % n, k = (int)ri(0, n - 1);
			P cand[] = {P(ri(-M, M), ri(-M, M)), h[i], (h[i] + h[j]) / 2, (h[i] + h[k]) / 2,
				h[i] + P(ri(-1, 1), ri(-1, 1)), h[i] * 2 - h[j]};
			for (P p : cand) if (abs(p.x) <= M && abs(p.y) <= M) checkHull(h, p);
		}
	}
	// hulls with many vertices
	rep(it,0,30) {
		vector<P> ps;
		int N = (int)ri(3, 3000);
		rep(i,0,N) { ll x = ri(-20000, 20000); ps.emplace_back(x, ri(0, 1) ? x * x : 2 * 20000LL * 20000 - x * x); }
		vector<P> h = convexHull(ps);
		rotate(h.begin(), h.begin() + ri(0, sz(h) - 1), h.end());
		int n = sz(h);
		rep(q,0,3000) {
			int i = (int)ri(0, n - 1), j = (i + 1) % n, k = (int)ri(0, n - 1);
			P cand[] = {P(ri(-20001, 20001), ri(-10, 2 * 20000LL * 20000 + 10)), h[i], (h[i] + h[j]) / 2,
				(h[i] + h[k]) / 2, h[i] + P(ri(-1, 1), ri(-1, 1))};
			for (P p : cand) checkHull(h, p);
		}
	}
}

int main() {
	testAgainstBrute();
	rep(it,0,100000) {
		int N = rand() % 15;
		vector<P> ps;
		rep(i,0,N) ps.emplace_back(rand() % 20 - 10, rand() % 20 - 10);
		vector<P> ps2 = convexHull(ps);
		if (ps2.empty()) continue;
		rep(it2,0,20) {
			int x = rand() % 22 - 11;
			int y = rand() % 22 - 11;
			P p{x,y};
			assert(inPolygon(ps2, p, true) == (inHull(ps2, p, true)));
			assert(inPolygon(ps2, p, false) == (inHull(ps2, p, false)));
		}
	}
	cout<<"Tests passed!"<<endl;
}
