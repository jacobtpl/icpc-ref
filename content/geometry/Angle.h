/**
 * Author: Simon Lindholm
 * Date: 2015-01-31
 * License: CC0
 * Source: me
 * Description: A class for ordering angles (as represented by integer points and
 *  a number of rotations around the origin). Useful for rotational sweeping.
 *  Sometimes also represents points or vectors.
 *  AngleT<T> takes int, long long or __int128; Angle is AngleT<ll>.
 *  Comparison is exact for all int and long long coordinates.
 *  +, - need the sums to fit. angleDiff needs |x|,|y| <= 2e9 for ll and
 *  returns an unreduced vector with coordinates up to 2 max(|x|,|y|)^2, so
 *  differences can be compared but not fed to angleDiff again above ~4e4.
 * Usage:
 *  vector<Angle> v = {w[0], w[0].t360() ...}; // sorted
 *  int j = 0; rep(i,0,n) { while (v[j] < v[i].t180()) ++j; }
 *  // sweeps j such that (j-i-1) represents the number of positively oriented triangles with vertices at 0 and i
 * Status: Used, works well
 */
#pragma once

#include "Point.h"

template<class T> struct AngleT {
	typedef AngleT Angle;
	T x, y;
	int t;
	AngleT(T x, T y, int t=0) : x(x), y(y), t(t) {}
	Angle operator-(Angle b) const { return {x-b.x, y-b.y, t}; }
	int half() const { assert(x || y); return y < 0 || (y == 0 && x < 0); }
	Angle t90() const { return {-y, x, t + (half() && x >= 0)}; }
	Angle t180() const { return {-x, -y, t + half()}; }
	Angle t360() const { return {x, y, t + 1}; }
};
typedef AngleT<ll> Angle;
#define AT template<class T, class Angle = AngleT<T>>
AT bool operator<(AngleT<T> a, AngleT<T> b) {
	// add a.dist2() and b.dist2() to also compare distances
	typedef typename Geo<T>::W W;
	return make_tuple(a.t, a.half(), W(a.y) * W(b.x)) <
	       make_tuple(b.t, b.half(), W(a.x) * W(b.y));
}
// Given two points, this calculates the smallest angle between
// them, i.e., the angle that covers the defined line segment.
AT pair<Angle, Angle> segmentAngles(AngleT<T> a, AngleT<T> b) {
	if (b < a) swap(a, b);
	return (b < a.t180() ?
	        make_pair(a, b) : make_pair(b, a.t360()));
}
AT Angle operator+(AngleT<T> a, AngleT<T> b) { // point a + vector b
	Angle r(a.x + b.x, a.y + b.y, a.t);
	if (a.t180() < r) r.t--;
	return r.t180() < a ? r.t360() : r;
}
AT Angle angleDiff(AngleT<T> a, AngleT<T> b) { // angle b - angle a
	int tu = b.t - a.t; a.t = b.t;
	return {a.x*b.x + a.y*b.y, a.x*b.y - a.y*b.x, tu - (b < a)};
}
#undef AT
