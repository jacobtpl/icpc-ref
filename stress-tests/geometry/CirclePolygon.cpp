#include "../utilities/template.h"

#include "../../content/geometry/CirclePolyIntersection.h"
typedef Point<double> P;
#include "../utilities/genPolygon.h"

// The reference loses half its digits near tangency, so it needs more
// than double: long double where that is wider, else __float128.
namespace orig{
#if LDBL_MANT_DIG >= 64
typedef long double R;
#else
typedef __float128 R;
R fabs(R x) { return gfabs(x); }
R sqrt(R x) { return gsqrt(x); }
R sin(R x) { return gsin(x); }
R atan2(R y, R x) { return gatan2(y, x); }
R acos(R x) { return atan2(sqrt((1 - x) * (1 + x)), x); }
R asin(R x) { return atan2(x, sqrt((1 - x) * (1 + x))); }
#endif
typedef Point<R> P;
R areaCT(P pa, P pb, R r) {
    if (pa.dist() < pb.dist()) swap(pa, pb);
    if (sgn(pb.dist()) == 0) return 0;
    R a = pb.dist(), b = pa.dist(), c = (pb - pa).dist();
    R sinB = fabs(pb.cross(pb - pa) / a / c), cosB = pb.dot(pb - pa) / a / c,
           sinC = fabs(pa.cross(pb) / a / b), cosC = pa.dot(pb) / a / b;
    R B = atan2(sinB, cosB), C = atan2(sinC, cosC);
    if (a > r) {
        R S = C / 2 * r * r, h = a * b * sinC / c;
        if (h < r && B < M_PI / 2)
            S -= (acos(h / r) * r * r - h * sqrt(r * r - h * h));
        return S;
    } else if (b > r) {
        R theta = M_PI - B - asin(sinB / r * a);
        return a * r * sin(theta) / 2 + (C - theta) / 2 * r * r;
    } else return sinC * a * b / 2;
}
R circlePoly(P c, R r, vector<P> poly) {
    R area = 0;
    rep(i,0,sz(poly)){
        auto a = poly[i] - c, b = poly[(i+1)%sz(poly)] - c;
        area += areaCT(a, b, r) * sgn(a.cross(b));
    }
    return area;
}
}

signed main() {
    double nearCentre = circlePoly(P(0,0), 1,
        {P(1e6,0), P(1e-9,1e-9), P(0,1e6)});
    assert(abs(nearCentre + acos(-1.) / 4) < 1e-8);
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

        if (abs(res1 - (double)res2) > 1e-8) {
            cout<<abs(res1-(double)res2)<<' '<<res1<<' '<<(double)res2<<endl;
            assert(false);
        }
    }
    cout<<"Tests passed!"<<endl;
}
