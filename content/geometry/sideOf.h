/**
 * Author: Ulf Lundstrom
 * Date: 2009-03-21
 * License: CC0
 * Source:
 * Description: Returns where p is as seen from s towards e.
 *  1/0/-1 <=> left/on line/right. The cross product is evaluated in
 *  the coordinate type T (watch out for overflow with int or
 *  long long) and its sign is taken with sgn: exact, unless
 *  geoEps<T> is set for a floating T (see Point.h).
 *  If the optional argument eps is given, 0 is returned if p is within
 *  distance eps from the line. Unlike geoEps this eps is a distance,
 *  in the units of the input. It has type F = Geo<T>::F (T itself for
 *  floating T), so no precision is lost with long double / __float128.
 * Usage:
 * 	typedef Point<double> P;
 * 	bool left = sideOf(p1,p2,q)==1;
 * 	bool near = sideOf(p1,p2,q,1e-9)==0; // within 1e-9 of the line
 * Status: tested
 */
#pragma once

#include "Point.h"

template<class P>
int sideOf(P s, P e, P p) { return sgn(s.cross(e, p)); }

template<class P>
int sideOf(const P& s, const P& e, const P& p, typename P::F eps) {
	typedef typename P::F F;
	F a = F((e-s).cross(p-s)), l = (e-s).dist()*eps;
	return (a > l) - (a < -l);
}
