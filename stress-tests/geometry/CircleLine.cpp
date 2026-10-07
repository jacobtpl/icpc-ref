#include "../utilities/template.h"
#include "../utilities/randGeo.h"

#include "../../content/geometry/lineDistance.h"
#include "../../content/geometry/CircleLine.h"

typedef Point<double> P;
ll irand(int lim) { return rand() % (2 * lim + 1) - lim; }

int main() {
    {
        auto res = circleLine(P(0, 0), 1, P(-1, -1), P(1, 1));
        assert(res.size() == 2);
        assert((res[1]-P(sqrt(2)/2, sqrt(2)/2)).dist() < 1e-8);
    }
    {
        auto res = circleLine(P(0, 0), 1, P(-5,  1), P(5, 1));
        assert(res.size() == 1);
        assert((res[0]-P(0,1)).dist() < 1e-8);
    }
    {
        auto res = circleLine(P(4, 4), 1, P(0,  0), P(5, 0));
        assert(res.size() == 0);
    }
	rep(it,0,100000) {
		P a = randIntPt(5);
		P b = randIntPt(5);
		P c = randIntPt(5);
		if (a == b) {
			// Not a well defined line
			continue;
		}
		double r = sqrt(rand() % 49);
		vector<P> points = circleLine(c, r, a, b);

		// Soundness
		assert(sz(points) <= 2);
		for (P p : points) {
			// Point is on circle
			assert(abs((p - c).dist() - r) < 1e-6);
			// Point is on line
			assert(lineDist(a, b, p) < 1e-6);
		}

		// Best-effort completeness check:
		// in some easy cases we must have points in the intersection.
		if ((a - c).dist() < r - 1e-6 || (b - c).dist() < r - 1e-6 || ((a + b) / 2 - c).dist() < r - 1e-6) {
			assert(!points.empty());
		}
	}
	// Exact number of intersections for integer input: compare the squared
	// distance s^2/|ab|^2 from the centre to the line with r^2.
	for (ll K : {1, 1000, 1000000}) rep(it,0,200000) {
		ll ax = irand(5) * K, ay = irand(5) * K, bx = irand(5) * K, by = irand(5) * K;
		ll cx = irand(5) * K, cy = irand(5) * K, ri = (rand() % 9) * K;
		if (ax == bx && ay == by) continue;
		P a((double)ax, (double)ay), b((double)bx, (double)by), c((double)cx, (double)cy);
		double r = (double)ri;
		__int128 s = (__int128)(bx - ax) * (cy - ay) - (__int128)(by - ay) * (cx - ax);
		__int128 d2 = (__int128)(bx - ax) * (bx - ax) + (__int128)(by - ay) * (by - ay);
		__int128 cmp = s * s - (__int128)ri * ri * d2;
		int expected = cmp > 0 ? 0 : cmp == 0 ? 1 : 2;
		vector<P> points = circleLine(c, r, a, b);
		assert(sz(points) == expected);
		for (P p : points) {
			assert(abs((p - c).dist() - r) < 1e-7 * (double)K);
			assert(lineDist(a, b, p) < 1e-7 * (double)K);
		}
		// Two results are distinct and ordered along the direction a -> b.
		if (expected == 2) assert((points[1] - points[0]).dot(b - a) > 0);
		// The answer only depends on the line, not on the two points chosen.
		vector<P> rev = circleLine(c, r, b, a);
		assert(sz(rev) == expected);
		rep(i,0,expected) assert((rev[i] - points[expected - 1 - i]).dist() < 1e-7 * (double)K);
	}
	// Non-integer input: every point returned is on both curves, and the line
	// is hit whenever its distance to the centre is clearly below r.
	rep(it,0,200000) {
		auto rd = []() { return (rand() % 20001 - 10000) / 100.0; };
		P a(rd(), rd()), b(rd(), rd()), c(rd(), rd());
		double r = abs(rd());
		if ((a - b).dist() < 1e-3) continue;
		vector<P> points = circleLine(c, r, a, b);
		double d = abs(lineDist(a, b, c));
		if (d < r - 1e-6) assert(sz(points) == 2);
		if (d > r + 1e-6) assert(sz(points) == 0);
		for (P p : points) {
			assert(abs((p - c).dist() - r) < 1e-6);
			assert(abs(lineDist(a, b, p)) < 1e-6);
		}
	}
    cout<<"Tests passed!"<<endl;
}
