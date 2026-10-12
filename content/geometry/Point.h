/**
 * Author: Ulf Lundstrom
 * Date: 2009-02-26
 * License: CC0
 * Source: My head with inspiration from tinyKACTL
 * Description: Class to handle points in the plane.
 * 	T can be e.g. double, long double or long long. (Avoid int.)
 * 	Integer T: dist2, dot and cross are evaluated in T, so with
 * 	long long keep $|x|,|y| \le 10^9$ ($2 \cdot 10^9$ for dist2/dist
 * 	of a single point); operator/ truncates. Lengths, angles, unit,
 * 	normal and rotate use F = T if T is floating, else double,
 * 	and never truncate.
 * 	Floating dist() is safe for all magnitudes.
 * 	mulDiv(a,b,c) is a*b/c with the product of integers kept in
 * 	\_\_int128 (truncates toward zero).
 * Status: Works fine, used a lot
 */
#pragma once

template <class T> int sgn(T x) { return (x > 0) - (x < 0); }
template<class T> T mulDiv(T a, T b, T c) {
	if constexpr (is_integral<T>::value)
		return T((__int128)a * b / c);
	else return a * b / c;
}
template<class T>
struct Point {
	typedef Point P;
	typedef conditional_t<is_floating_point_v<T>,T,double> F;
	typedef Point<F> PF;
	T x, y;
	explicit Point(T x=0, T y=0) : x(x), y(y) {}
	bool operator<(P p) const { return tie(x,y) < tie(p.x,p.y); }
	bool operator==(P p) const { return tie(x,y)==tie(p.x,p.y); }
	bool operator!=(P p) const { return !(*this == p); }
	P& operator+=(P p) { x += p.x, y += p.y; return *this; }
	P& operator-=(P p) { x -= p.x, y -= p.y; return *this; }
	P operator+(P p) const { return P(x+p.x, y+p.y); }
	P operator-(P p) const { return P(x-p.x, y-p.y); }
	P operator*(T d) const { return P(x*d, y*d); }
	P operator/(T d) const { return P(x/d, y/d); }
	T dot(P p) const { return x*p.x + y*p.y; }
	T cross(P p) const { return x*p.y - y*p.x; }
	T cross(P a, P b) const { return (a-*this).cross(b-*this); }
	T dist2() const { return x*x + y*y; }
	F dist() const { // hypot only if dist2() over/underflows F
		F d = F(dist2()), e = numeric_limits<F>::epsilon();
		if (!is_floating_point_v<T> || isnormal(d*e) || (!x && !y))
			return sqrt(d);
		return hypot(F(x), F(y));
	}
	// angle to x-axis in interval [-pi, pi]
	F angle() const { return atan2(F(y), F(x)); }
	PF unit() const { // makes dist()=1
		F d = dist(); return PF(F(x)/d, F(y)/d); }
	P perp() const { return P(-y, x); } // rotates +90 degrees
	PF normal() const { return perp().unit(); }
	// returns point rotated 'a' radians ccw around the origin
	PF rotate(F a) const {
		F c = cos(a), s = sin(a), u = F(x), v = F(y);
		return PF(u*c - v*s, u*s + v*c); }
	friend ostream& operator<<(ostream& os, P p) {
		return os << "(" << p.x << "," << p.y << ")"; }
};
