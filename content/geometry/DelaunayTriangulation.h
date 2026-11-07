/**
 * Author: Mattias de Zalenski
 * Date: Unknown
 * Source: Geometry in C
 * Description: Computes the Delaunay triangulation of a set of points.
 *  Each circumcircle contains none of the input points.
 *  Calls trifun(i,j,k) with the indices of every triangle, counter-clockwise.
 *  If any three points are collinear or any four are on the same circle, behavior is undefined
 *  (use FastDelaunay for such input).
 *  Integer coordinates are lifted exactly in \_\_int128: valid for absolute values up to $10^9$.
 *  Floating point coordinates are centred and lifted in long double, which is not exact.
 * Time: O(n^2)
 * Status: stress-tested
 */
#pragma once

#include "Point.h"
#include "3dHull.h"

template<class P, class F>
void delaunay(vector<P>& ps, F trifun) {
	typedef decltype(ps[0].x) T;
	typedef conditional_t<is_integral<T>::value,
		__int128, long double> S;
	if (sz(ps) == 3) { int d = (ps[0].cross(ps[1], ps[2]) < 0);
		trifun(0,1+d,2-d); }
	if (sz(ps) <= 3) return;
	P lo = ps[0], hi = lo; // centre: keeps |coords| and helps doubles
	for (P p : ps) lo = P(min(lo.x, p.x), min(lo.y, p.y)),
		hi = P(max(hi.x, p.x), max(hi.y, p.y));
	P o = (lo + hi) / 2;
	vector<Point3D<S>> p3;
	for (P p : ps) { S x = p.x - o.x, y = p.y - o.y;
		p3.emplace_back(x, y, x*x + y*y); }
	for (auto t : hull3d(p3)) if ((p3[t.b]-p3[t.a]).
			cross(p3[t.c]-p3[t.a]).z < 0)
		trifun(t.a, t.c, t.b);
}
