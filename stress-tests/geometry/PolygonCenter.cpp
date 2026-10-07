#include "../utilities/template.h"

#include "../../content/geometry/PolygonCenter.h"

typedef long double ld;
typedef __int128 lll;
mt19937_64 rng(11);
ll ri(ll a, ll b) { return uniform_int_distribution<ll>(a, b)(rng); }

// Oracle: fan triangulation from v[0] in exact integer arithmetic. The
// centroid is (cx / (3*A2), cy / (3*A2)) where A2 is twice the signed area.
struct Exact { lll cx, cy, A2; };
Exact oracle(const vector<Point<ll>>& v) {
	Exact r{0, 0, 0};
	rep(i,1,sz(v)-1) {
		lll a = (lll)(v[i].x-v[0].x) * (v[i+1].y-v[0].y) - (lll)(v[i].y-v[0].y) * (v[i+1].x-v[0].x);
		r.A2 += a;
		r.cx += a * (v[0].x + v[i].x + v[i+1].x);
		r.cy += a * (v[0].y + v[i].y + v[i+1].y);
	}
	return r;
}

int main() {
	// known shapes
	assert((polygonCenter({P(0,0), P(4,0), P(4,2), P(0,2)}) - P(2,1)).dist() < 1e-12);
	assert((polygonCenter({P(0,0), P(3,0), P(0,3)}) - P(1,1)).dist() < 1e-12);
	// L-shape: 2x2 square minus the unit square in the top right corner
	assert((polygonCenter({P(0,0), P(2,0), P(2,1), P(1,1), P(1,2), P(0,2)})
		- P(5/6.0, 5/6.0)).dist() < 1e-12);

	ll tested = 0;
	rep(it,0,1500000) {
		int n = (int)ri(3, 8);
		ll lim = it % 3 == 0 ? 1000000 : it % 3 == 1 ? ri(1, 5) : 1000;
		ll ox = it % 4 == 0 ? ri(-lim, lim) * 3 : 0, oy = it % 4 == 0 ? ri(-lim, lim) * 3 : 0;
		vector<Point<ll>> q(n);
		for (auto& p : q) p = Point<ll>(ox + ri(-lim, lim), oy + ri(-lim, lim));
		if (it % 5 == 0) q[ri(0, n-1)] = q[ri(0, n-1)]; // repeated vertex
		Exact e = oracle(q);
		if (e.A2 == 0) continue; // no area, centre undefined
		tested++;
		vector<P> v;
		for (auto p : q) v.emplace_back((double)p.x, (double)p.y);
		P c = polygonCenter(v);
		ld wx = (ld)e.cx / (3 * (ld)e.A2), wy = (ld)e.cy / (3 * (ld)e.A2);
		// cancellation: absolute error is about eps * scale^3 / area
		ld sc = (ld)(4 * lim);
		ld tol = 1e-9L * max<ld>(1, sc) * max<ld>(1, sc * sc / fabsl((ld)e.A2));
		if (!(fabsl(c.x - wx) <= tol && fabsl(c.y - wy) <= tol)) {
			cout << setprecision(17) << "got " << c << " want (" << (double)wx << "," << (double)wy << ")" << endl;
			for (auto p : v) cout << p << ' ';
			cout << endl;
			abort();
		}
		// orientation and starting vertex must not matter
		vector<P> w = v;
		reverse(all(w));
		rotate(w.begin(), w.begin() + ri(0, n-1), w.end());
		P c2 = polygonCenter(w);
		assert(abs(c.x - c2.x) <= 2*tol && abs(c.y - c2.y) <= 2*tol);
	}
	assert(tested > 1000000);
	cout<<"Tests passed!"<<endl;
}
