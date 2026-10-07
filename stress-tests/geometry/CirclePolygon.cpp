#include "../utilities/template.h"

#include "../../content/geometry/CirclePolyIntersection.h"
#include "../utilities/genPolygon.h"

namespace orig{
typedef Point<long double> P;
long double areaCT(P pa, P pb, long double r) {
    if (pa.dist() < pb.dist()) swap(pa, pb);
    if (sgn(pb.dist()) == 0) return 0;
    long double a = pb.dist(), b = pa.dist(), c = (pb - pa).dist();
    long double sinB = fabs(pb.cross(pb - pa) / a / c), cosB = pb.dot(pb - pa) / a / c,
           sinC = fabs(pa.cross(pb) / a / b), cosC = pa.dot(pb) / a / b;
    long double B = atan2(sinB, cosB), C = atan2(sinC, cosC);
    if (a > r) {
        long double S = C / 2 * r * r, h = a * b * sinC / c;
        if (h < r && B < M_PI / 2)
            S -= (acos(h / r) * r * r - h * sqrt(r * r - h * h));
        return S;
    } else if (b > r) {
        long double theta = M_PI - B - asin(sinB / r * a);
        return a * r * sin(theta) / 2 + (C - theta) / 2 * r * r;
    } else return sinC * a * b / 2;
}
long double circlePoly(P c, long double r, vector<P> poly) {
    long double area = 0;
    rep(i,0,sz(poly)){
        auto a = poly[i] - c, b = poly[(i+1)%sz(poly)] - c;
        area += areaCT(a, b, r) * sgn(a.cross(b));
    }
    return area;
}
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);
    const int lim=5;
    for (int i=0;i<100000; i++) {

        vector<Point<int>> pts;
        for (int j=0; j<10; j++) {
            int x = rand()%lim, y = rand()%lim;
            pts.push_back(Point<int>(x, y));
        }

        auto polyInt = genPolygon(pts);

        int cx = rand()%lim, cy = rand()%lim;
        auto c = P(cx, cy);
        auto c2 = orig::P(cx, cy);
        double r= rand()%(2*lim);

        vector<P> poly;
        vector<orig::P> poly2;
        for (auto j: polyInt) {
            poly.push_back(P(j.x, j.y));
            poly2.push_back(orig::P(j.x, j.y));
        }
        auto res1 = circlePoly(c, r, poly);
        auto res2 = orig::circlePoly(c2, r, poly2);

        if (abs(res1 - res2) > 1e-8) {
            cout<<abs(res1-res2)<<' '<<res1<<' '<<res2<<endl;
            assert(false);
        }
    }
	// Closed forms and invariances.
	const double pi = acos(-1);
	assert(circlePoly(P(1, 2), 3, {}) == 0);
	for (int i=0;i<100000; i++) {
		auto rd = [](double lim) { return (rand() % 20001 - 10000) / 10000.0 * lim; };
		P c(rd(100), rd(100));
		double r = abs(rd(10)) + 1e-3, L = 10 + abs(rd(50));
		// circle strictly inside a square
		vector<P> sq = {c + P(-L, -L), c + P(L, -L), c + P(L, L), c + P(-L, L)};
		assert(abs(circlePoly(c, r, sq) - pi * r * r) < 1e-9 * r * r);
		// clockwise polygons give the negated area
		assert(abs(circlePoly(c, r, {sq[3], sq[2], sq[1], sq[0]}) + pi * r * r) < 1e-9 * r * r);
		// centre on an edge / at a corner of the square: half / quarter disc
		vector<P> half = {c + P(0, -L), c + P(L, -L), c + P(L, L), c + P(0, L)};
		assert(abs(circlePoly(c, r, half) - pi * r * r / 2) < 1e-9 * r * r);
		vector<P> quarter = {c, c + P(L, 0), c + P(L, L), c + P(0, L)};
		assert(abs(circlePoly(c, r, quarter) - pi * r * r / 4) < 1e-9 * r * r);
		// circle far away from the polygon, and radius 0
		assert(abs(circlePoly(c + P(3 * L, 0), r, sq)) < 1e-9);
		assert(abs(circlePoly(c, 0, sq)) < 1e-9);
		// polygon inside the circle: plain polygon area
		vector<Point<int>> pts;
		for (int j=0; j<8; j++) pts.push_back(Point<int>(rand()%21-10, rand()%21-10));
		auto polyInt = genPolygon(pts);
		vector<P> poly, moved, scaled;
		double K = 1 + rand() % 1000;
		P off(rd(1000), rd(1000)), c0(rd(10), rd(10));
		for (auto p : polyInt) {
			poly.push_back(P(p.x, p.y));
			moved.push_back(P(p.x, p.y) + off);
			scaled.push_back(P(p.x, p.y) * K);
		}
		double area = polygonArea2(poly) / 2;
		assert(abs(circlePoly(c0, 100, poly) - area) < 1e-7);
		// translating / scaling everything translates / scales the answer
		double r0 = abs(rd(15)), res = circlePoly(c0, r0, poly);
		assert(-1e-9 <= res && res <= min(area, pi * r0 * r0) + 1e-9);
		assert(abs(circlePoly(c0 + off, r0, moved) - res) < 1e-7);
		assert(abs(circlePoly(c0 * K, r0 * K, scaled) - res * K * K) < 1e-9 * K * K);
	}
    cout<<"Tests passed!"<<endl;
}
