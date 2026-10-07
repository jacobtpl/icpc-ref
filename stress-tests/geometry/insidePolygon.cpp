#include <bits/stdc++.h>
using namespace std;

#define rep(i, a, b) for(int i = a; i < int(b); ++i)
#define all(x) x.begin(), x.end()
#define sz(x) (int)(x).size()

typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;


const double EPS =1e-8;
#include "../utilities/genPolygon.h"
#include "../../content/geometry/InsidePolygon.h"
namespace old {

#include "../../content/geometry/OnSegment.h"
#include "../../content/geometry/SegmentDistance.h"

template<class It, class P>
bool insidePolygon(It begin, It end, const P& p,
		bool strict = true) {
	int n = 0; //number of isects with line from p to (inf,p.y)
	for (It i = begin, j = end-1; i != end; j = i++) {
		//if p is on edge of polygon
		if (onSegment(*i, *j, p)) return !strict;
		//or: if (segDist(*i, *j, p) <= epsilon) return !strict;
		//increment n if segment intersects line from p
		n += (max(i->y,j->y) > p.y && min(i->y,j->y) <= p.y &&
				((*j-*i).cross(p-*i) > 0) == (i->y <= p.y));
	}
	return n&1; //inside if odd number of intersections
}
}
typedef Point<double> P;
bool eq(P a, P b) {
    return (a-b).dist()<EPS;
}
const int NUMPOLY=100;
const int PTPERPOLY=100;
void test(int numPts, int range) {
    rep(i,0,NUMPOLY) {
        vector<P> poly;
        rep(j,0, numPts)
            poly.push_back(P(rand()%range, rand()%range));
        poly = genPolygon(poly);
        rep(i,0,PTPERPOLY){
            P p(rand()%range, rand()%range);
            assert(inPolygon(poly, p, true) == old::insidePolygon(all(poly), p, true));
            assert(inPolygon(poly, p, false) == old::insidePolygon(all(poly), p, false));
        }
    }

}
// Independent oracle on integer points: exact boundary test, otherwise the
// winding number from the sum of signed angles.
// Returns 0 = outside, 1 = on the boundary, 2 = strictly inside.
typedef Point<ll> PL;
typedef __int128_t lll;
int oracle(const vector<PL>& poly, PL a) {
	int n = sz(poly); long double ang = 0;
	rep(i,0,n) {
		PL u = poly[i] - a, v = poly[(i+1)%n] - a;
		lll cr = (lll)u.x*v.y - (lll)u.y*v.x, dt = (lll)u.x*v.x + (lll)u.y*v.y;
		if (cr == 0 && dt <= 0) return 1;
		ang += atan2l((long double)cr, (long double)dt);
	}
	ang = fabsl(ang);
	assert(ang < 0.1 || fabsl(ang - 2*acosl(-1)) < 0.1); // polygon is simple
	return ang > 3 ? 2 : 0;
}
void check(vector<PL> poly, PL a) {
	int o = oracle(poly, a);
	assert(inPolygon(poly, a, true) == (o == 2));
	assert(inPolygon(poly, a, false) == (o >= 1));
}
// all rotations and both orientations, against a list of query points
void checkAll(vector<PL> poly, const vector<PL>& qs, bool allRot = true) {
	rep(rev,0,2) {
		rep(r,0,allRot ? max(sz(poly), 1) : 2) {
			for (PL a : qs) check(poly, a);
			if (sz(poly)) rotate(poly.begin(), poly.begin() + 1, poly.end());
		}
		reverse(all(poly));
	}
}
ll sgnl(lll x) { return (x > 0) - (x < 0); }
bool segsTouch(PL a, PL b, PL c, PL d) { // closed segments intersect
	auto cr = [](PL o, PL p, PL q) { return sgnl((lll)(p.x-o.x)*(q.y-o.y) - (lll)(p.y-o.y)*(q.x-o.x)); };
	auto on = [&](PL s, PL e, PL p) { return cr(p, s, e) == 0 && (s-p).dot(e-p) <= 0; };
	if (on(a,b,c) || on(a,b,d) || on(c,d,a) || on(c,d,b)) return true;
	return cr(a,b,c) * cr(a,b,d) < 0 && cr(c,d,a) * cr(c,d,b) < 0;
}
bool isSimple(const vector<PL>& p) {
	int n = sz(p);
	rep(i,0,n) rep(j,0,i) {
		int i2 = (i+1)%n, j2 = (j+1)%n;
		if (i2 == j || j2 == i) { // adjacent edges only share their common vertex
			PL s = i2 == j ? p[i] : p[j], m = i2 == j ? p[j] : p[i], e = i2 == j ? p[j2] : p[i2];
			if ((lll)(s-m).cross(e-m) == 0 && (s-m).dot(e-m) > 0) return false;
		} else if (segsTouch(p[i], p[i2], p[j], p[j2])) return false;
	}
	return true;
}
// star-shaped around (0,0): one point on each of n distinct directions
vector<PL> starPolygon(int n, ll C) {
	set<pair<ll,ll>> dirs;
	while (sz(dirs) < n) {
		ll x = rand() % 21 - 10, y = rand() % 21 - 10, g = __gcd(abs(x), abs(y));
		if (g) dirs.insert({x / g, y / g});
	}
	vector<PL> v;
	for (auto d : dirs) {
		ll k = rand() % (C / 10) + 1;
		v.push_back(PL(d.first * k, d.second * k));
	}
	auto half = [](PL p) { return p.y < 0 || (p.y == 0 && p.x < 0); };
	sort(all(v), [&](PL a, PL b) {
		return half(a) != half(b) ? half(a) < half(b) : a.cross(b) > 0; });
	return v;
}
// rectilinear histogram with collinear vertices and horizontal edges
vector<PL> histogram(int w, int H) {
	vector<PL> v = {PL(0,0), PL(w,0)};
	for (int x = w; x > 0; x--) {
		ll h = rand() % H + 1;
		v.push_back(PL(x, h)); v.push_back(PL(x-1, h));
	}
	return v;
}
void testInteger() {
	// degenerate polygons
	vector<PL> grid;
	rep(x,-3,4) rep(y,-3,4) grid.push_back(PL(x,y));
	checkAll({}, grid);
	checkAll({PL(1,-1)}, grid);
	checkAll({PL(-2,-1), PL(2,1)}, grid);
	checkAll({PL(0,0), PL(2,0), PL(2,2), PL(0,2)}, grid);
	checkAll({PL(-2,0), PL(0,-2), PL(2,0), PL(0,2)}, grid);
	// the usage example from the header
	{
		vector<PL> v = {PL(4,4), PL(1,2), PL(2,1)};
		assert(inPolygon(v, PL(3,3), false));
		assert(inPolygon(v, PL(3,3), true));
		assert(inPolygon(v, PL(4,4), false) && !inPolygon(v, PL(4,4), true));
		assert(!inPolygon(v, PL(1,1), false));
		checkAll(v, grid);
	}
	// small random simple polygons, every grid point as query
	int simple = 0;
	rep(it,0,40000) {
		int n = rand() % 9 + 1, C = rand() % 5 + 2;
		vector<Point<double>> pd;
		rep(i,0,n) pd.push_back(Point<double>(rand() % C, rand() % C));
		pd = genPolygon(pd);
		vector<PL> poly;
		for (auto p : pd) poly.push_back(PL((ll)p.x, (ll)p.y));
		if (!isSimple(poly)) continue;
		simple++;
		vector<PL> qs;
		rep(x,-1,C+1) rep(y,-1,C+1) qs.push_back(PL(x,y));
		checkAll(poly, qs, it % 8 == 0);
	}
	assert(simple > 20000);
	rep(it,0,3000) {
		vector<PL> poly = it % 2 ? histogram(rand() % 8 + 1, rand() % 4 + 1)
			: starPolygon(rand() % 10 + 3, 60);
		vector<PL> qs;
		rep(x,-7,10) rep(y,-7,8) qs.push_back(PL(x,y));
		checkAll(poly, qs, it % 8 == 0);
	}
	// coordinates up to 1e9 (products up to 8e18 fit in ll)
	rep(it,0,3000) {
		const ll C = 1000000000;
		vector<PL> poly = starPolygon(rand() % 30 + 3, C);
		if (it % 3 == 0) { // rectangle touching the limits
			ll a = C - rand() % 3, b = C - rand() % 3;
			poly = {PL(-a,-b), PL(a,-b), PL(a,b), PL(-a,b)};
		}
		vector<PL> qs;
		int n = sz(poly);
		rep(i,0,n) {
			PL p = poly[i], q = poly[(i+1)%n];
			qs.push_back(p);
			rep(dx,-1,2) rep(dy,-1,2) qs.push_back(p + PL(dx,dy));
			qs.push_back(PL((p.x + q.x) / 2, (p.y + q.y) / 2));
			qs.push_back(PL(rand() % (2*C+1) - C, p.y)); // ray through a vertex
			qs.push_back(PL(rand() % (2*C+1) - C, rand() % (2*C+1) - C));
			qs.push_back(PL(-C, p.y)); qs.push_back(PL(C, p.y));
		}
		for (auto& a : qs) a.x = max(-C, min(C, a.x)), a.y = max(-C, min(C, a.y));
		checkAll(poly, qs, false);
	}
}

int main() {
    srand(7);
    testInteger();
    test(20,5);
    test(1001,100);
    test(1000,1000);
    cout<<"Tests passed!"<<endl;
}
