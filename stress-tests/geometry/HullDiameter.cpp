#include "../utilities/template.h"

#include "../../content/geometry/ConvexHull.h"
#include "../../content/geometry/HullDiameter.h"

ll brute(const vector<P>& ps) {
	ll r = 0;
	rep(i,0,sz(ps)) rep(j,0,i) r = max(r, (ps[i] - ps[j]).dist2());
	return r;
}
// checks every rotation of the hull (any starting vertex is valid input)
void check(vector<P> hull, bool allRot = true) {
	if (hull.empty()) return;
	ll want = brute(hull);
	rep(r,0,allRot ? sz(hull) : 3) {
		auto pa = hullDiameter(hull);
		assert((pa[0] - pa[1]).dist2() == want);
		assert(count(all(hull), pa[0]) && count(all(hull), pa[1]));
		rotate(hull.begin(), hull.begin() + (allRot ? 1 : rand() % sz(hull)), hull.end());
	}
}

int main() {
	srand(2);
	rep(it,0,1000000) {
		int N = (rand() % 10) + 1;
		vector<Point<ll>> ps;
		rep(i,0,N) {
			ps.emplace_back(rand() % 11 - 5, rand() % 11 - 5);
		}
		ll r1 = 0;
		rep(i,0,N) rep(j,0,i) {
			r1 = max(r1, (ps[i] - ps[j]).dist2());
		}
		auto pa = hullDiameter(convexHull(ps));
		ll r2 = ps.empty() ? 0LL : (pa[0] - pa[1]).dist2();
		assert(r1 == r2);
	}
	// all rotations, tiny inputs
	rep(it,0,300000) {
		int N = rand() % 12 + 1, C = rand() % 6 + 1;
		vector<P> ps;
		rep(i,0,N) ps.emplace_back(rand() % (2*C+1) - C, rand() % (2*C+1) - C);
		check(convexHull(ps));
	}
	// explicit shapes
	check({P(7,-3)});
	check({P(0,0), P(5,0)});
	check({P(0,0), P(1,0), P(0,1)});
	check({P(0,0), P(1,0), P(1,1), P(0,1)}); // two equal diagonals
	check({P(0,0), P(1000000000,0), P(1000000000,1), P(0,1)});
	check({P(-1000000000,-1000000000), P(1000000000,-1000000000),
			P(1000000000,1000000000), P(-1000000000,1000000000)});
	// larger hulls: random, on a circle, on a parabola, thin, near 1e9
	rep(it,0,3000) {
		int N = rand() % 300 + 1, type = it % 5;
		ll C = it % 2 ? 1000000000 : 1000;
		vector<P> ps;
		rep(i,0,N) {
			ll x = rand() % (2*C+1) - C, y = rand() % (2*C+1) - C;
			if (type == 1) {
				double a = rand(), R = (double)C;
				x = llround(R * cos(a)), y = llround(R * sin(a));
			}
			if (type == 2) x = x % 30000, y = x * x - C;
			if (type == 3) y %= 3;
			if (type == 4) x = x / 1000 * 1000, y = y / 1000 * 1000;
			ps.emplace_back(x, y);
		}
		check(convexHull(ps), it % 10 == 0);
	}
	cout<<"Tests passed!"<<endl;
}
