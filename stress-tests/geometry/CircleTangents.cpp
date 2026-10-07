#include "../utilities/template.h"

#include "../../content/geometry/CircleTangents.h"
#include "../../content/geometry/lineDistance.h"
#include "../utilities/randGeo.h"

typedef Point<double> P;

ll irand(int lim) { return rand() % (2 * lim + 1) - lim; }

signed main() {
    for (int i = 0; i < 1000000; i++) {
        P c1 = randIntPt(5), c2 = randIntPt(5);
        double r1 = sqrt(rand()%20), r2 = sqrt(rand()%20);
        for (auto sgn : {-1, 1}) {
            auto tans = tangents(c1, r1, c2, sgn * r2);

            if (tans.size() ==1) {
                assert((tans[0].first - tans[0].second).dist() < 1e-8);
                assert(abs((tans[0].first-c1).dist() - r1) < 1e-8);
                assert(abs((tans[0].first-c2).dist() - r2) < 1e-8);
            } else if (tans.size() == 2) {
                for (auto l : tans) {
                    assert(abs(abs(lineDist(l.first, l.second, c1))-r1) < 1e-8);
                    assert(abs(abs(lineDist(l.first, l.second, c2))-r2) < 1e-8);
                }
            }
        }
    }
	// Exact number of tangents for integer input, and the geometry of each.
	for (ll K : {1, 1000, 1000000}) rep(it,0,300000) {
		ll x1 = irand(5) * K, y1 = irand(5) * K, x2 = irand(5) * K, y2 = irand(5) * K;
		ll a1 = (rand() % 8) * K, a2 = (rand() % 8) * K;
		P c1((double)x1, (double)y1), c2((double)x2, (double)y2);
		double r1 = (double)a1, r2 = (double)a2, eps = 1e-7 * (double)K;
		ll d2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
		for (int sgn : {-1, 1}) {
			// external: 2 unless one circle is inside the other (1 if they touch
			// internally); internal: 2 unless the circles overlap (1 if they
			// touch externally); none for concentric circles.
			ll dr = a1 - sgn * a2, h2 = d2 - dr * dr;
			int expected = d2 == 0 || h2 < 0 ? 0 : h2 == 0 ? 1 : 2;
			auto tans = tangents(c1, r1, c2, sgn * r2);
			assert(sz(tans) == expected);
			for (auto l : tans) {
				P u = l.first - c1, v = l.second - c2, t = l.second - l.first;
				assert(abs(u.dist() - r1) < eps);
				assert(abs(v.dist() - r2) < eps);
				// the tangent is perpendicular to both radii
				assert(abs(u.dot(t)) < eps * (double)K * 20);
				assert(abs(v.dot(t)) < eps * (double)K * 20);
				// radii are parallel (external) or anti-parallel (internal)
				assert(abs(u.cross(v)) < eps * (double)K * 20);
				if (a1 && a2) assert((u.dot(v) > 0) == (sgn == 1));
			}
			if (expected == 1) assert((tans[0].first - tans[0].second).dist() < eps);
			if (expected == 2) {
				if (a1) assert((tans[0].first - tans[1].first).dist() > eps);
				if (a2) assert((tans[0].second - tans[1].second).dist() > eps);
			}
			// Tangents from a point: r2 = 0 gives the point itself as .second.
			if (a2 == 0) for (auto l : tans) assert(l.second == c2);
		}
	}
    cout<<"Tests passed!"<<endl;
}
