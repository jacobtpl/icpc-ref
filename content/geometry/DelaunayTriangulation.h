/**
 * Author: Mattias de Zalenski
 * Date: Unknown
 * Source: Geometry in C
 * Description: Computes the Delaunay triangulation of a set of points.
 *  Each circumcircle contains none of the input points.
 *  Calls trifun(i,j,k) with the indices of every triangle, counter-clockwise.
 *  If any three points are collinear or any four are on the same circle, behavior is undefined
 *  (use FastDelaunay for such input).
 *  Coordinates are centred and lifted in W = Geo<T>::W (see Point.h).
 *  Integer T is exact while 4th powers of the coordinate range fit in W:
 *  Point<ll> (W = __int128) is valid for absolute values up to 1e9,
 *  Point<int> (W = ll) up to 1e4.
 *  Floating T is not exact; by default it is lifted in T itself.
 * Usage:
 *  vector<Point<ll>> ps = ...;
 *  delaunay(ps, [&](int i, int j, int k) { tris.push_back({i,j,k}); });
 * Time: O(n^2)
 * Status: stress-tested
 */
#pragma once

#include "Point.h"
#include "3dHull.h"

template<class P, class F>
void delaunay(vector<P>& ps, F trifun) {
	typedef decltype(ps[0].x) T;
	typedef typename Geo<T>::W S;
	if (sz(ps) == 3) { int d = (ps[0].cross(ps[1], ps[2]) < 0);
		trifun(0,1+d,2-d); }
	if (sz(ps) <= 3) return;
	P lo = ps[0], hi = lo; // centre: keeps |coords| and helps doubles
	for (P p : ps) lo = P(min(lo.x, p.x), min(lo.y, p.y)),
		hi = P(max(hi.x, p.x), max(hi.y, p.y));
	P o = (lo + hi) / 2;
	vector<Point3D<S>> p3;
	for (P p : ps) { S x = S(p.x) - S(o.x), y = S(p.y) - S(o.y);
		p3.emplace_back(x, y, x*x + y*y); }
	for (auto t : hull3d(p3)) if ((p3[t.b]-p3[t.a]).
			cross(p3[t.c]-p3[t.a]).z < 0)
		trifun(t.a, t.c, t.b);
}
