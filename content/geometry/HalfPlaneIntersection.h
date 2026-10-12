/**
 * Author: Benq
 * Description: Returns vertices of half-plane intersection.
   * A half-plane is the area to the left of a ray, which is defined
   * by a point \texttt{p} and a direction \texttt{dp}. 
   * Area of intersection should be sufficiently precise when all inputs
   * are integers with magnitude $\le 10^5$. Intersection must be bounded.
   * Probably works with floating point too (but EPS might need to be adjusted?).
   * Use RayT<Fraction> for exact rational intersections with zero tolerance.
 * Time: O(N\log N)
 * Source: Own
 * Status: tested on
	* https://open.kattis.com/problems/bigbrother 
	* https://icpc.kattis.com/problems/domes
	* https://codeforces.com/gym/100962 B
 */

#include "Point.h"

using P = Point<double>;
using vP = vector<P>;
const double EPS = 1e-9; // adjust as needed

// -1 if lower half, 0 if origin, 1 if upper half, needs to be int type
template<class P>
int half(P x) { return x.y != 0 ? sgn(x.y) : -sgn(x.x); }
template<class P>
bool angleCmp(P a, P b) { int A = half(a), B = half(b);
    return A == B ? a.cross(b) > 0 : A < B; }

template<class T> struct RayT {
	using P = Point<T>;
	P p, dp; // origin, direction
	P isect(const RayT& L) const {
		return p+dp*(L.dp.cross(L.p-p)/L.dp.cross(dp)); }
	bool operator<(const RayT& L) const {
		return angleCmp(dp,L.dp); }
};
using Ray = RayT<double>;

template<class T>
vector<Point<T>> halfPlaneIsect(vector<RayT<T>> rays, bool add_bounds = false) {
	using P = Point<T>;
	using Ray = RayT<T>;
	T eps = 0;
	if constexpr (is_floating_point<T>::value) eps = T(EPS);
	if (add_bounds) { // bound input by rectangle [0,DX] x [0,DY]
		int DX = 1e9, DY = 1e9;
		rays.pb({P{0,0},P{1,0}});
		rays.pb({P{DX,0},P{0,1}});
		rays.pb({P{DX,DY},P{-1,0}});
		rays.pb({P{0,DY},P{0,-1}});
	}
	sort(all(rays)); // sort rays by angle
	{ // remove parallel rays
		vector<Ray> nrays;
		for(auto& t: rays) {
			if (!sz(nrays) || nrays.back().dp.cross(t.dp) > eps) { nrays.pb(t); continue; }
			// last two rays are parallel, keep only one
			if (t.dp.cross(t.p-nrays.back().p) > 0) nrays.back() = t;
		}
		swap(rays, nrays);
	}
	auto bad = [&](const Ray& a, const Ray& b, const Ray& c) {
		P p1 = a.isect(b), p2 = b.isect(c);
		if ((p2-p1).dot(b.dp) <= eps) { /// tolerance for coincident intersections
			if (a.dp.cross(c.dp) <= 0) return 2; // isect(a,b,c) = empty
			return 1; // isect(a,c) == isect(a,b,c)
		}
		return 0; // all three rays matter
	};
	#define reduce(t) \
		while (sz(poly) > 1) { \
			int b = bad(poly.at(sz(poly)-2),poly.back(),t); \
			if (b == 2) return {}; \
			if (b == 1) poly.pop_back(); \
			else break; \
		}
	deque<Ray> poly;
    for(auto& t: rays) { reduce(t); poly.pb(t); }
	for(;;poly.pop_front()) {
		reduce(poly[0]);
		if (!bad(poly.back(),poly[0],poly[1])) break;
	}
	assert(sz(poly) >= 3); // expect nonzero area
	vector<P> poly_points; rep(i,0,sz(poly))
		poly_points.pb(poly[i].isect(poly[(i+1)%sz(poly)]));
	return poly_points;
}
