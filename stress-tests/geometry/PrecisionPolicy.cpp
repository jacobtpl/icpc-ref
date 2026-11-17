// Every geometry header that computes non-integer values must work for
// double, long double and __float128 (without -lquadmath) and for int,
// long long and __int128, in both -std=c++20 and -std=gnu++20, with
// accuracy that follows the declared type. Also checks the geoEps hook.
#include "../utilities/template.h"

#include "../../content/geometry/Point.h"
#include "../../content/geometry/Point3D.h"
#include "../../content/geometry/CircleLine.h"
#include "../../content/geometry/CircleTangents.h"
#include "../../content/geometry/CircleIntersection.h"
#include "../../content/geometry/CirclePolyIntersection.h"
#include "../../content/geometry/MinimumEnclosingCircle.h"
#include "../../content/geometry/sphericalDistance.h"
#include "../../content/geometry/sideOf.h"
#include "../../content/geometry/HalfPlaneIntersection.h"
#include "../../content/geometry/PolygonUnion.h"
#include "../../content/geometry/PolygonCenter.h"
#include "../../content/geometry/PolygonCut.h"
#include "../../content/geometry/ClosestPair.h"
#include "../../content/geometry/FastDelaunay.h"
#include "../../content/geometry/ManhattanMST.h"
#include "../../content/geometry/kdTree.h"
#include "../../content/geometry/Angle.h"
#include "../../content/geometry/lineIntersection.h"
#include "../../content/geometry/SegmentDistance.h"
#include "../../content/geometry/ConvexHull.h"

typedef __float128 Q;
typedef long double ld;
template<class T> ld err(T a, T b) { return (ld)gfabs(a - b); }

// u = unit roundoff of T; results must be within tol*u
template<class T> void floats(ld u) {
	typedef Point<T> P;
	static_assert(is_same_v<typename P::F, T>);
	static_assert(is_same_v<decltype(P().dist()), T>);
	T two = 2, s2 = gsqrt(two), pi = gatan2(T(0), T(-1));
	assert(err(s2 * s2, two) < 8 * u);
	assert(err(P(3, 4).dist(), T(5)) == 0);
	assert(err(P(1, 1).angle() * 4, pi) < 8 * u);
	assert(err(P(1, 0).rotate(pi / 2).y, T(1)) < 8 * u);
	assert(err(Point3D<T>(1, 2, 2).dist(), T(3)) == 0);
	assert(err(Point3D<T>(0, 0, 1).rotate(pi, Point3D<T>(1, 0, 0)).z, T(-1)) < 8 * u);
	// the radius keeps the precision of T: r = 1 + 2^-k is a secant
	T r = 1 + T(u) * 4;
	assert(r > 1 && sz(circleLine(P(0, 0), r, P(-1, 1), P(1, 1))) == 2);
	assert(sz(circleLine(P(0, 0), T(1), P(-1, 1), P(1, 1))) == 1);
	auto cl = circleLine(P(0, 0), T(5), P(-10, 3), P(10, 3));
	assert(sz(cl) == 2 && err(cl[0].x, T(-4)) < 40 * u && err(cl[1].x, T(4)) < 40 * u);
	assert(sz(tangents(P(0, 0), T(1), P(4, 0), T(1))) == 2);
	assert(sz(tangents(P(0, 0), T(1), P(2, 0), T(-1))) == 1);
	pair<P, P> io;
	assert(circleInter(P(0, 0), P(4, 0), T(3), T(3), &io));
	assert(err(io.first.x, T(2)) < 40 * u && err(io.second.y * io.second.y, T(5)) < 80 * u);
	assert(!circleInter(P(0, 0), P(4, 0), T(1), T(2), &io));
	vector<P> sq = {P(0, 0), P(2, 0), P(2, 2), P(0, 2)};
	assert(err(circlePoly(P(0, 0), T(1), sq), pi / 4) < 40 * u);
	assert(err(ccRadius(P(0, 0), P(2, 0), P(0, 2)), s2) < 40 * u);
	// thin triangle: compensated determinants keep the relative error
	P cc = ccCenter(P(0, 0), P(T(1e6), 1), P(T(2e6), 3));
	assert(err((cc - P(0, 0)).dist() / (cc - P(T(2e6), 3)).dist(), T(1)) < 100 * u);
	auto m = mec(vector<P>{P(0, 0), P(2, 0), P(1, 1)}, T(u * 1e4));
	assert(err(m.second, T(1)) < 1e5 * u && err(m.first.x, T(1)) < 1e5 * u);
	assert(err(sphericalDistance(T(0), pi / 2, pi / 2, pi / 2, T(1)), pi / 2) < 40 * u);
	assert(sideOf(P(0, 0), P(1, 0), P(0, T(u))) == 1);
	assert(sideOf(P(0, 0), P(1, 0), P(0, T(u)), T(u * 2)) == 0);
	vector<RayT<T>> rs = {{P(0,0),P(1,0)}, {P(1,0),P(0,1)}, {P(0,1),P(-1,0)}, {P(0,1),P(0,-1)}};
	auto hp = halfPlaneIsect(rs, false, T(u * 512));
	static_assert(is_same_v<decltype(hp), vector<P>>);
	assert(sz(hp) == 4);
	assert(sz(halfPlaneIsect(vector<RayT<T>>{}, true, T(u * 512), 7)) == 4);
	vector<vector<P>> un = {sq, {P(1, 1), P(3, 1), P(3, 3), P(1, 3)}};
	assert(err(polyUnion(un), T(7)) < 40 * u);
	assert(err(polygonCenter(sq).x, T(1)) < 8 * u);
	assert(sz(polygonCut(sq, P(1, 0), P(1, 2))) == 4);
	assert(err(segDist(P(0, 0), P(2, 0), P(1, 1)), T(1)) == 0);
	assert(lineInter(P(0, 0), P(2, 2), P(0, 2), P(2, 0)).second == P(1, 1));
	assert(closest(sq).first == P(0, 0));
	{ // fractional and duplicate coordinates
		vector<P> v, dup = {P(0, T(1e16)), P(0, T(1e16))};
		rep(i,0,300) v.push_back(P(T(i * 37 % 300) / 512, T(i * 91 % 300) / 512));
		T best = 9;
		rep(i,0,300) rep(k,0,i) best = min(best, (v[i] - v[k]).dist2());
		auto c = closest(v);
		assert((c.first - c.second).dist2() == best);
		assert(closest(dup).first == dup[0]);
		// large offset: p - d rounds, brute force must still agree
		rep(it,0,200) {
			vector<P> w;
			rep(i,0,40) w.push_back(P(T(1e15) + T((i * 7919 + it * 31) % 97) / 4,
				T(-3e15) + T((i * 104729 + it * 17) % 89) / 4));
			best = T(1e30);
			rep(i,0,40) rep(k,0,i) best = min(best, (w[i] - w[k]).dist2());
			c = closest(w);
			assert((c.first - c.second).dist2() == best);
		}
	}
	assert(sz(triangulate(sq)) == 6);
	assert(sz(convexHull(sq)) == 4);
	assert(KDTree<P>(sq).nearest(P(3, 3)).first == 2);
	assert(sz(manhattanMST(sq)) >= 3);
}

template<class T, class F, class W> void ints(T big) {
	typedef Point<T> P;
	static_assert(is_same_v<typename P::F, F>);
	static_assert(is_same_v<typename Geo<T>::W, W>);
	static_assert(geoInt<T> && !geoFloat<T>);
	assert(P(3, 4).dist() == F(5));
	assert(Point3D<T>(1, 2, 2).dist() == F(3));
	assert(sgn(T(-7)) == -1 && sgn(T(0)) == 0);
	// product in W, not T
	assert(mulDiv(big, big, big) == big);
	AngleT<T> a(big, 1), b(big, 2);
	assert(a < b && !(b < a));
	vector<P> sq = {P(0, 0), P(big, 0), P(big, big), P(0, big)};
	// polygonCenter(Point<int>) used to overflow and truncate
	auto c = polygonCenter(sq);
	static_assert(is_same_v<decltype(c), Point<F>>);
	assert(c.x == F(big) / 2 && c.y == F(big) / 2);
	auto cut = polygonCut(vector<P>{P(0, 0), P(1, 0), P(1, 1), P(0, 1)}, P(0, 0), P(2, 1));
	static_assert(is_same_v<decltype(cut), vector<Point<F>>>);
	assert(sz(cut) == 3 && count(all(cut), Point<F>(1, F(1) / 2)) == 1);
	vector<P> s2 = {P(0, 0), P(2, 0), P(2, 2), P(0, 2)};
	vector<vector<P>> un = {s2, {P(1, 1), P(3, 1), P(3, 3), P(1, 3)}};
	assert(polyUnion(un) == F(7)); // was 28 for __int128 in -std=c++20
	vector<RayT<T>> rs = {{P(0,0),P(1,0)}, {P(2,0),P(0,1)}, {P(0,2),P(-1,0)}, {P(0,2),P(0,-1)}};
	auto hp = halfPlaneIsect(rs);
	static_assert(is_same_v<decltype(hp), vector<Point<F>>>);
	assert(sz(hp) == 4);
	rep(i,0,4) assert(hp[i].x == F(i == 1 || i == 2 ? 2 : 0));
	assert((halfPlaneIsect<T, ld>(rs)).size() == 4);
	{ // the vertex is computed in the requested type, not in F
		vector<RayT<T>> tr = {{P(0,0),P(3,1)}, {P(3,1),P(-3,6)}, {P(0,7),P(0,-7)}};
		auto g = halfPlaneIsect<T, ld>(tr);
		assert(sz(g) == 3 && g[1].y == 1 && g[1].x == 3);
		vector<RayT<T>> th = {{P(0,0),P(1,0)}, {P(1,0),P(-1,3)}, {P(0,3),P(0,-1)},
			{P(0,1),P(1,0)}};
		for (auto p : halfPlaneIsect<T, ld>(th)) if (p.y == 1 && p.x > 0)
			assert(fabsl(p.x - 2 / 3.0L) < 1e-18L);
	}
	assert(closest(s2).first == P(0, 0));
	assert(sz(triangulate(s2)) == 6);
	assert(KDTree<P>(s2).nearest(P(3, 3)).first == 2);
	assert(sz(manhattanMST(s2)) >= 3);
	assert(lineInter(P(0, 0), P(2, 2), P(0, 2), P(2, 0)).second == P(1, 1));
	assert(sideOf(P(0, 0), P(1, 0), P(0, 1)) == 1);
}

int main() {
	static_assert(geoFloat<Q> && geoInt<__int128>); // also in strict mode
	floats<double>(0x1p-53L);
	floats<ld>(LDBL_EPSILON / 2);
	floats<Q>(0x1p-113L);
	// __float128 really carries ~34 digits
	Q s = gsqrt(Q(2));
	assert(s * s != Q(2) && gfabs(s * s - 2) < Q(1e-33));
	assert((ld)s == sqrtl(2.0L) && Q((ld)s) != s);

	ints<int, double, ll>(40000);
	ints<ll, double, __int128>(3000000000LL);
	ints<__int128, ld, __int128>(3000000000LL);

	// float: lengths whose square is subnormal or overflows
	assert(fabsf(Point<float>(3e-23f, 0).dist() / 3e-23f - 1) < 1e-6f);
	assert(fabsf(Point3D<float>(0, 3e-23f, 0).dist() / 3e-23f - 1) < 1e-6f);
	assert(fabsf(Point<float>(3e30f, 4e30f).dist() / 5e30f - 1) < 1e-6f);
	assert(fabs(Point<double>(3e-160, 4e-160).dist() / 5e-160 - 1) < 1e-15);

	// geoEps: off by default, then a tolerance for floating sgn only
	typedef Point<double> P;
	assert(sgn(1e-300) == 1 && sideOf(P(0, 0), P(1, 0), P(0.5, 1e-12)) == 1);
	geoEps<double> = 1e-9;
	assert(sgn(1e-10) == 0 && sgn(-1e-10) == 0 && sgn(1e-8) == 1 && sgn(-1e-8) == -1);
	assert(sideOf(P(0, 0), P(1, 0), P(0.5, 1e-12)) == 0);
	assert(sz(circleLine(P(0, 0), 1.0, P(-1, 1 + 1e-12), P(1, 1 + 1e-12))) == 1);
	assert(sgn(1) == 1 && sgn(1e-10L) == 1); // other types unaffected
	geoEps<double> = 0;
	assert(sgn(1e-10) == 1);
	cout << "Tests passed!" << endl;
}
