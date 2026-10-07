#include "../utilities/template.h"

#include "../../content/geometry/DelaunayTriangulation.h"
#define ll double
#include "../../content/geometry/ConvexHull.h"
#undef ll
#include "../../content/geometry/PolygonArea.h"
#include "../../content/geometry/circumcircle.h"

typedef Point<double> P;
int main() {
	feenableexcept(29);
	rep(it,0,100000) {{
		vector<P> ps;
		int N = rand() % 20 + 1;
		rep(i,0,N) {
			ps.emplace_back(rand() % 100 - 50, rand() % 100 - 50);
		}

		auto coc = [&](int i, int j, int k, int l) {
			double a = (ps[i] - ps[j]).dist();
			double b = (ps[j] - ps[k]).dist();
			double c = (ps[k] - ps[l]).dist();
			double d = (ps[l] - ps[i]).dist();
			double e = (ps[i] - ps[k]).dist();
			double f = (ps[j] - ps[l]).dist();
			double q = a*c + b*d - e*f;
			return abs(q) < 1e-4;
		};

		rep(i,0,N) rep(j,0,i) rep(k,0,j) {
			if (ps[i].cross(ps[j], ps[k]) == 0) {  goto fail; }
		}
		rep(i,0,N) rep(j,0,i) rep(k,0,j) rep(l,0,k) {
			if (coc(i,j,k,l) || coc(i,j,l,k) || coc(i,l,j,k) || coc(i,l,k,j)) { goto fail; }
		}

		auto fail = [&]() {
			cout << "Points:" << endl;
			for(auto &p: ps) {
				cout << p.x << ' ' << p.y << endl;
			}

			cout << "Triangles:" << endl;
			delaunay(ps, [&](int i, int j, int k) {
				cout << i << ' ' << j << ' ' << k << endl;
			});

			abort();
		};

		double sumar = 0;
		vi used(N);
		delaunay(ps, [&](int i, int j, int k) {
			used[i] = used[j] = used[k] = 1;
			double ar = ps[i].cross(ps[j], ps[k]);
			if (ar < -1e-4) fail();
			sumar += ar;
			P c = ccCenter(ps[i], ps[j], ps[k]);
			double ra = ccRadius(ps[i], ps[j], ps[k]);
			rep(l,0,N) {
				if ((ps[l] - c).dist() < ra - 1e-5) fail();
			}
		});
		if (N >= 3) rep(i,0,N) if (!used[i]) fail();

		vector<P> hull = convexHull(ps);
		double ar2 = polygonArea2(hull);
		if (abs(sumar - ar2) > 1e-4) fail();

		continue; }
fail:;
	}
	// Fewer than three points: no triangles; exactly three: one ccw triangle.
	rep(n,0,3) {
		vector<P> ps(n);
		rep(i,0,n) ps[i] = P(i, i * i);
		delaunay(ps, [&](int, int, int) { assert(false); });
	}
	// Exact check on integer points in general position, up to |x| <= 1e6:
	// ccw triangles with empty circumcircles covering the hull, every point used.
	typedef __int128 lll;
	auto inCircle = [](Point<ll> a, Point<ll> b, Point<ll> c, Point<ll> d) {
		a = a - d, b = b - d, c = c - d;
		return (lll)a.dist2() * b.cross(c) + (lll)b.dist2() * c.cross(a) + (lll)c.dist2() * a.cross(b);
	};
	for (int K : {3, 10, 1000, 1000000}) rep(it,0,20000) {
		int N = rand() % 12 + 3;
		if (K == 3) N = rand() % 3 + 3;
		vector<Point<ll>> q;
		if (it % 2) rep(i,0,N) q.emplace_back(rand() % (2*K+1) - K, rand() % (2*K+1) - K);
		else { // thin triangles
			Point<ll> d(rand() % K + 1, rand() % K + 1);
			int m = max(1, K / (int)max(d.x, d.y));
			rep(i,0,N) q.push_back(d * (rand() % (2*m+1) - m) + Point<ll>(rand() % 5 - 2, rand() % 5 - 2));
		}
		bool bad = 0;
		rep(i,0,N) rep(j,0,i) rep(k,0,j) {
			if (q[i].cross(q[j], q[k]) == 0) bad = 1;
			else rep(l,0,k) if (inCircle(q[i], q[j], q[k], q[l]) == 0) bad = 1;
		}
		if (bad) continue;
		vector<P> ps;
		for (auto p : q) ps.emplace_back((double)p.x, (double)p.y);
		lll sum = 0, area = 0;
		vi used(N);
		delaunay(ps, [&](int i, int j, int k) {
			used[i] = used[j] = used[k] = 1;
			ll ar = q[i].cross(q[j], q[k]);
			assert(ar > 0);
			sum += ar;
			rep(l,0,N) assert(inCircle(q[i], q[j], q[k], q[l]) <= 0);
		});
		rep(i,0,N) assert(used[i]);
		vector<P> hull = convexHull(ps);
		rep(i,0,sz(hull)) area += (lll)llround(hull[i].cross(hull[(i + 1) % sz(hull)]));
		assert(sum == area);
	}
	cout<<"Tests passed!"<<endl;
}
