/**
 * Author: Ulf Lundstrom with inspiration from tinyKACTL
 * Date: 2009-04-14
 * License: CC0
 * Source:
 * Description: Class to handle points in 3D space.
 * 	T can be double, long double, __float128, int, long long, __int128
 * 	or Fraction; see Point.h. dot, cross and dist2 are evaluated in T
 * 	(long long: keep coordinates <= 1e9). Lengths, angles, unit,
 * 	normal and rotate return a floating type and never truncate.
 * Usage:
 * 	typedef Point3D<double> P3;
 * 	P3 a(1, 0, 0), b(0, 1, 0);
 * 	P3 n = a.cross(b);            // (0,0,1)
 * 	double len = (a - b).dist();  // sqrt(2)
 * 	auto r = a.rotate(acos(-1) / 2, n); // (0,1,0), type Point3D<F>
 * Status: tested, except for phi and theta
 */
#pragma once

#include "Point.h"

template<class T> struct Point3D {
	typedef Point3D P;
	typedef const P& R;
	typedef typename Geo<T>::F F;
	typedef Point3D<F> PF;
	T x, y, z;
	explicit Point3D(T x=0, T y=0, T z=0) : x(x), y(y), z(z) {}
	bool operator<(R p) const {
		return tie(x, y, z) < tie(p.x, p.y, p.z); }
	bool operator==(R p) const {
		return tie(x, y, z) == tie(p.x, p.y, p.z); }
	bool operator!=(R p) const { return !(*this == p); }
	P operator+(R p) const { return P(x+p.x, y+p.y, z+p.z); }
	P operator-(R p) const { return P(x-p.x, y-p.y, z-p.z); }
	P operator*(T d) const { return P(x*d, y*d, z*d); }
	P operator/(T d) const { return P(x/d, y/d, z/d); }
	T dot(R p) const { return x*p.x + y*p.y + z*p.z; }
	P cross(R p) const {
		return P(y*p.z - z*p.y, z*p.x - x*p.z, x*p.y - y*p.x);
	}
	T dist2() const { return x*x + y*y + z*z; }
	F dist() const { // hypot only if dist2() over/underflows F
		F d = F(dist2()), big = F(sizeof(F) < 8 ? 1e30 : 1e290);
		if (!geoFloat<T> || (d > 1 / big && d < big) ||
			(x == 0 && y == 0 && z == 0)) return gsqrt(d);
		F a = gfabs(F(x)), b = gfabs(F(y)), c = gfabs(F(z)),
		  m = max(a, max(b, c));
		a /= m, b /= m, c /= m;
		return m * gsqrt(a*a + b*b + c*c);
	}
	//Azimuthal angle (longitude) to x-axis in interval [-pi, pi]
	F phi() const { return gatan2(F(y), F(x)); }
	//Zenith angle (latitude) to the z-axis in interval [0, pi]
	F theta() const { return gatan2(ghypot(F(x), F(y)), F(z)); }
	PF unit() const { //makes dist()=1
		F d = dist(); return PF(F(x)/d, F(y)/d, F(z)/d); }
	//returns unit vector normal to *this and p
	PF normal(P p) const { return cross(p).unit(); }
	//returns point rotated 'angle' radians ccw around axis
	PF rotate(F angle, P axis) const {
		F s = gsin(angle), c = gcos(angle);
		PF u = axis.unit(), v = PF(F(x), F(y), F(z));
		return u*v.dot(u)*(1-c) + v*c - v.cross(u)*s;
	}
};
