#include "../utilities/template.h"

#include "../../content/geometry/PolygonArea.h"
#include "../../content/geometry/PolygonCut.h"
#include "../../content/geometry/sideOf.h"
#include "../../content/geometry/InsidePolygon.h"
#include "../../content/geometry/SegmentIntersection.h"

typedef Point<double> P;
int main() {
	vector<P> thin = {P(1e7+0.3,1e7+0.3+1e-6),
		P(1e7+10.7,1e7+10.7-1e-6), P(1e7+11.4,1e7+11.4+1e-6)};
	auto clipped = polygonCut(thin, P(1e7-1.7,1e7-1.7), P(1e7+12.3,1e7+12.3));
	assert(sz(clipped) == 3);
	assert((clipped[0] - P(1e7+5.5,1e7+5.5)).dist() < 1e-7);
	assert(clipped[1] == thin[1]);
	assert((clipped[2] - P(1e7+11.05,1e7+11.05)).dist() < 1e-7);
	rep(it,0,500) {
		int N = rand() % 10 + 3;
		vector<P> ps;
		rep(i,0,N) ps.emplace_back(rand() % 10 - 5, rand() % 10 - 5);
		P p(rand() % 10 - 5), q(rand() % 10 - 5);
		rep(i,0,N) rep(j,i+1,N) {
			P a = ps[i], b = ps[(i+1)%N];
			P c = ps[j], d = ps[(j+1)%N];
			P r1, r2;
			auto r = segInter(a, b, c, d);
			if (sz(r) == 2) goto fail;
			if (sz(r) == 1) {
				if (i+1 == j || (j+1) % N == i) ;
				else goto fail;
			}
		}
		if (p == q) { fail: continue; }

		int count = 0;
		const int ITS = 400000;
		rep(it,0,ITS) {
			double x = rand() / (RAND_MAX + 1.0) * 10 - 5;
			double y = rand() / (RAND_MAX + 1.0) * 10 - 5;
			if (!inPolygon(ps, P{x,y}, true)) continue;
			if (sideOf(p, q, P{x,y}) > 0) continue;
			count++;
		}
		double approxArea = (double)count / ITS * 100;

		ps = polygonCut(ps, p, q);
		double realArea = ps.empty() ? 0.0 : abs(polygonArea2(ps) / 2.0);

		// cout << setprecision(2) << fixed;
		assert(realArea - approxArea < 2e-1);
		// cout << N << ' ' << realArea << '\t' << approxArea << '\t' << realArea - approxArea << endl;

		// cerr << N << endl;
		// for(auto &x: ps) {
			// cout << x.x << ' ' << x.y << endl;
		// }
	}
	cout<<"Tests passed!"<<endl;
}
