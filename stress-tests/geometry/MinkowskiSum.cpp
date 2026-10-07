#include "../utilities/template.h"
#define pb push_back

#include "../../content/geometry/ConvexHull.h"
#include "../../content/geometry/MinkowskiSum.h"

mt19937 rng(12345);
int ri(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }

// convex hull in CCW order that keeps collinear boundary points
vector<P> hullCol(vector<P> p) {
	sort(all(p)); p.erase(unique(all(p)), p.end());
	if (sz(p) < 3) return p;
	vector<P> lo, up;
	for (P q : p) {
		while (sz(lo) >= 2 && lo[sz(lo)-2].cross(lo.back(), q) < 0) lo.pop_back();
		lo.push_back(q);
	}
	if (sz(lo) == sz(p)) return convexHull(p); // all collinear
	for (int i = sz(p) - 1; i >= 0; i--) {
		while (sz(up) >= 2 && up[sz(up)-2].cross(up.back(), p[i]) < 0) up.pop_back();
		up.push_back(p[i]);
	}
	lo.pop_back(); up.pop_back();
	lo.insert(lo.end(), all(up));
	return lo;
}

void fail(const char* why, const vector<P>& a, const vector<P>& b, const vector<P>& got) {
	cerr << "FAIL (" << why << ")\na:";
	for (P x : a) cerr << ' ' << x;
	cerr << "\nb:";
	for (P x : b) cerr << ' ' << x;
	cerr << "\ngot:";
	rep(i,0,min(sz(got), 20)) cerr << ' ' << got[i];
	cerr << endl;
	exit(1);
}

// strictHull: both inputs are strictly convex with >= 3 vertices (or a point),
// then the output must be exactly the strictly convex hull of all pairwise sums.
void check(const vector<P>& a, const vector<P>& b, bool strict) {
	vector<P> s;
	for (P x : a) for (P y : b) s.push_back(x + y);
	vector<P> want = convexHull(s);
	vector<P> got = minkowski_sum(a, b);
	if (a.empty() || b.empty()) {
		if (!got.empty()) fail("empty", a, b, got);
		return;
	}
	if (sz(got) > sz(a) + sz(b)) fail("too many points", a, b, got);
	if (convexHull(got) != want) fail("wrong hull", a, b, got);
	// output is a convex CCW polygon without repeated vertices
	int n = sz(got);
	if (n > 1) rep(i,0,n) {
		if (got[i] == got[(i + 1) % n]) fail("repeated vertex", a, b, got);
		if (got[i].cross(got[(i + 1) % n], got[(i + 2) % n]) < 0) fail("not convex", a, b, got);
	}
	if (strict && sz(want) >= 3) {
		rotate(got.begin(), min_element(all(got)), got.end());
		if (got != want) fail("not the strict hull in CCW order", a, b, got);
	}
}

vector<P> randPts(int n, int C) {
	vector<P> p;
	rep(i,0,n) p.emplace_back(ri(-C, C), ri(-C, C));
	return p;
}
void rot(vector<P>& h) {
	if (sz(h)) rotate(h.begin(), h.begin() + ri(0, sz(h) - 1), h.end());
}

int main() {
	// strictly convex polygons, points and segments (as 2-gons), any start vertex
	rep(it,0,300000) {
		int C = ri(1, 6);
		vector<P> a = convexHull(randPts(ri(0, 8), C)), b = convexHull(randPts(ri(0, 8), C));
		rot(a), rot(b);
		check(a, b, sz(a) != 2 && sz(b) != 2);
	}
	// convex polygons with collinear points on the boundary
	rep(it,0,300000) {
		int C = ri(1, 4);
		vector<P> a = hullCol(randPts(ri(1, 10), C)), b = hullCol(randPts(ri(1, 10), C));
		rot(a), rot(b);
		check(a, b, false);
	}
	// segments (2-gons) against polygons with collinear boundary points
	rep(it,0,300000) {
		int C = ri(1, 4);
		vector<P> a = convexHull(randPts(2, C)), b = hullCol(randPts(ri(1, 12), C));
		rot(a), rot(b);
		check(a, b, false);
		check(b, a, false);
	}
	// a segment plus a polygon with collinear boundary points and a parallel edge
	{
		vector<P> a = {P(-1,2), P(1,2)};
		vector<P> b = {P(0,0), P(-2,0), P(-3,-2), P(-2,-3), P(2,-1), P(3,0)};
		check(a, b, false);
		check(b, a, false);
	}
	// large coordinates (|x|,|y| <= 1e9)
	rep(it,0,20000) {
		int C = 1000000000;
		vector<P> a = convexHull(randPts(ri(1, 12), C)), b = convexHull(randPts(ri(1, 12), C));
		rot(a), rot(b);
		check(a, b, sz(a) != 2 && sz(b) != 2);
	}
	// many vertices: points on a parabola-like convex curve
	rep(it,0,20) {
		auto gen = [&](int n) {
			vector<P> p;
			rep(i,0,n) { ll x = ri(-2000, 2000); p.emplace_back(x, ri(0, 1) ? x * x : 8000000 - x * x); }
			return convexHull(p);
		};
		vector<P> a = gen(ri(1, 300)), b = gen(ri(1, 300));
		rot(a), rot(b);
		check(a, b, sz(a) != 2 && sz(b) != 2);
	}
	cout<<"Tests passed!"<<endl;
}
