/**
 * Author: Ulf Lundstrom with inspiration from tinyKACTL
 * Date: 2009-04-14
 * License: CC0
 * Source:
 * Description: Class to handle points in 3D space.
 * 	T can be e.g. double, long double or long long.
 * 	Integer T: dot, cross and dist2 are evaluated in T (with long long
 * 	keep coordinates of the operands $\le 10^9$) and operator/ truncates.
 * 	Lengths, angles, unit, normal and rotate use F = T if T is
 * 	floating, else double, and never truncate.
 * Usage:
 * Status: tested, except for phi and theta
 */
#pragma once

template<class T> struct Point3D {
	typedef Point3D P;
	typedef const P& R;
	typedef conditional_t<is_floating_point_v<T>,T,double> F;
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
		F d = F(dist2()), e = numeric_limits<F>::epsilon();
		if (!is_floating_point_v<T> || isnormal(d*e) ||
			(!x && !y && !z)) return sqrt(d);
		return hypot(F(x), F(y), F(z));
	}
	//Azimuthal angle (longitude) to x-axis in interval [-pi, pi]
	F phi() const { return atan2(F(y), F(x)); }
	//Zenith angle (latitude) to the z-axis in interval [0, pi]
	F theta() const { return atan2(hypot(F(x), F(y)), F(z)); }
	PF unit() const { //makes dist()=1
		F d = dist(); return PF(F(x)/d, F(y)/d, F(z)/d); }
	//returns unit vector normal to *this and p
	PF normal(P p) const { return cross(p).unit(); }
	//returns point rotated 'angle' radians ccw around axis
	PF rotate(F angle, P axis) const {
		F s = sin(angle), c = cos(angle);
		PF u = axis.unit(), v = PF(F(x), F(y), F(z));
		return u*v.dot(u)*(1-c) + v*c - v.cross(u)*s;
	}
};
