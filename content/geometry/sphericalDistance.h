/**
 * Author: Ulf Lundstrom
 * Date: 2009-04-07
 * License: CC0
 * Source: My geometric reasoning
 * Description: Returns the shortest distance on the sphere with radius radius between the points with azimuthal angles (longitude) f1 (phi1) and f2 (phi2) from x axis and zenith angles (latitude) t1 (theta1) and t2 (theta2) from z axis (0 = north pole). All angles measured in radians. The algorithm starts by converting the spherical coordinates to cartesian coordinates so if that is what you have you can use only the two last rows. d and s are the lengths of the difference and the sum of the two unit vectors; the angle between them is 2 atan2(d, s), which stays accurate (and never NaN) for nearby and antipodal points alike.
 * Everything is computed in the type T of the first argument (double,
 * long double or __float128; an integer first argument gives double),
 * which is also the result type. Name it to be explicit:
 * sphericalDistance<long double>(...).
 * Usage:
 *  double pi = acos(-1.0);
 *  // equator, a quarter turn apart, unit sphere: pi/2
 *  double d = sphericalDistance(0.0, pi/2, pi/2, pi/2, 1.0);
 *  long double e = sphericalDistance<long double>(0, 1, 2, 1, 6371);
 * Status: tested on kattis:airlinehub
 */
#pragma once

#include "Point.h"

template<class A, class T = GeoR<A>>
T sphericalDistance(A f1, type_identity_t<T> t1,
		type_identity_t<T> f2, type_identity_t<T> t2,
		type_identity_t<T> radius) {
	T x1 = gsin(t1)*gcos(T(f1)), y1 = gsin(t1)*gsin(T(f1)), z1 = gcos(t1);
	T x2 = gsin(t2)*gcos(f2), y2 = gsin(t2)*gsin(f2), z2 = gcos(t2);
	T dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;
	T ax = x2 + x1, ay = y2 + y1, az = z2 + z1;
	T d = gsqrt(dx*dx + dy*dy + dz*dz);
	T s = gsqrt(ax*ax + ay*ay + az*az);
	return radius*2*gatan2(d, s);
}
