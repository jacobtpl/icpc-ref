/**
 * Author: Ulf Lundstrom
 * Date: 2009-02-26
 * License: CC0
 * Source: My head with inspiration from tinyKACTL
 * Description: Class to handle points in the plane. Pick precision
 *  by the point type; all geometry headers follow:
 *   typedef Point<double> P; // or long double, __float128 (no
 *                            // -lquadmath needed), int, ll, __int128
 *  dot, cross, dist2 are computed in T (ll: keep |x|,|y| <= 1e9).
 *  dist, angle, unit, rotate return a floating type and never truncate.
 *  Epsilon: exact by default. For a tolerance on floats set, in main,
 *   geoEps<double> = 1e-7;
 *  It is compared to cross/dot products: use ~1000 * u * M^2 for
 *  coordinates up to M (u = 2e-16 double, 1e-19 long double).
 *  Integer input: use an integer T and no eps.
 * Usage:
 *  typedef Point<double> P;
 *  P a(1, 2), b(4, 6);
 *  double d = (b - a).dist();      // 5
 *  auto c = a.cross(b);            // 1*6 - 2*4 = -2, type T
 *  Point<ll> q(3, 4); double l = q.dist(); // F of long long is double
 * Status: Works fine, used a lot
 */
#pragma once

// Number types per coordinate type T. Keep all rows; edit one only to
// change what Point<T> computes with:
//  F = type of lengths and non-integer results
//  W = type algorithms form products of dot/cross/dist2 results in
//      (dot, cross, dist2 themselves return T)
#define GEO(T, FT, WT) \
	template<> struct Geo<T> { typedef FT F; typedef WT W; };
template<class T> struct Geo { // any other type, e.g. Fraction
	typedef conditional_t<is_class_v<T>, double, T> F; typedef T W; };
#ifndef GEO_DOUBLE_W // long double: wider products for double points
#define GEO_DOUBLE_W double
#endif
//  T            F            W
GEO(double,      double,      GEO_DOUBLE_W)
GEO(long double, long double, long double)
GEO(__float128,  __float128,  __float128)
GEO(int,         double,      ll)
GEO(ll,          double,      __int128)
GEO(__int128,    long double, __int128) // nothing wider: stay < 2^63
#undef GEO
template<class T> T geoEps = T(0); // tolerance of sgn(); 0 = exact

template<class T> constexpr bool geoFloat =
	is_floating_point_v<T> || is_same_v<T, __float128>;
template<class T> constexpr bool geoInt =
	is_integral_v<T> || is_same_v<T, __int128>;
// type of non-integer results: F for integer T, else T itself
template<class T> using GeoR =
	conditional_t<geoInt<T>, typename Geo<T>::F, T>;
// type such results are computed in before being stored as GeoR<T>
template<class T> using GeoC =
	conditional_t<geoInt<T>, typename Geo<T>::F, typename Geo<T>::W>;
template<class T> int sgn(T x) {
	if constexpr (geoFloat<T>)
		return (x > geoEps<T>) - (x < -geoEps<T>);
	else return (x > 0) - (x < 0);
}
// __float128 math without libquadmath: glibc's functions on Linux,
// else (macOS; or define GEO_OWN_F128) plain C++, accurate to a few
// ulp for angles |x| < 1e12.
namespace q128 {
typedef __float128 Q;
#if defined(__GLIBC__) && !defined(GEO_OWN_F128)
#define GEO_Q(f) inline Q f(Q x) { return f##f128(x); }
GEO_Q(sqrt) GEO_Q(sin) GEO_Q(cos) GEO_Q(fabs)
#undef GEO_Q
inline Q atan2(Q y, Q x) { return atan2f128(y, x); }
inline Q hypot(Q x, Q y) { return hypotf128(x, y); }
inline Q fma(Q a, Q b, Q c) { return fmaf128(a, b, c); }
#else
inline Q fabs(Q x) { return x < 0 ? -x : x; }
inline Q sqrt(Q x) {
	if (!(x > 0) || x == x + x) return x < 0 ? (x - x) / (x - x) : x;
	Q s = 1; // scale into double range, then Newton from a double
	for (; x > 0x1p1000; x *= 0x1p-200) s *= 0x1p100;
	for (; x < 0x1p-1000; x *= 0x1p200) s *= 0x1p-100;
	Q r = std::sqrt((double)x);
	r = (r + x / r) / 2; r = (r + x / r) / 2;
	return r * s; }
inline Q hypot(Q x, Q y) {
	x = fabs(x), y = fabs(y); if (x < y) swap(x, y);
	if (x == 0 || x == x + x) return x;
	y /= x; return x * sqrt(1 + y * y); }
inline void sincos(Q x, Q& s, Q& c) {
	static const array<Q, 34> inv = [] { // 1/i!
		array<Q, 34> a; a[0] = 1;
		rep(i,1,34) a[i] = a[i-1] / i;
		return a; }();
	ll k = (ll)(x / 0x1.921fb54442d18p+0 + (x < 0 ? -.5 : .5));
	Q n = k, r = x - n * 0x1.921fb54442d18p+0 // x - k*pi/2
		- n * 0x1.1a62633145cp-54 - n * 0x1.b839a252049c1p-104,
		r2 = -r * r, S = 0, C = 0;
	for (int i = 32; i >= 0; i -= 2) // Taylor series, |r| <= pi/4
		C = C * r2 + inv[i], S = S * r2 + inv[i+1];
	S *= r; s = k & 1 ? C : S, c = k & 1 ? -S : C;
	if (k & 2) s = -s, c = -c; }
inline Q sin(Q x) { Q s, c; sincos(x, s, c); return s; }
inline Q cos(Q x) { Q s, c; sincos(x, s, c); return c; }
inline Q atan2(Q y, Q x) { // double atan2 plus one Newton step
	Q h = hypot(x, y), s, c;
	if (!(h > 0) || h == h + h) return std::atan2((double)y, (double)x);
	x /= h, y /= h; Q t = std::atan2((double)y, (double)x);
	sincos(t, s, c); Q d = y * c - x * s;
	return t + (d + d * d * d / 6); }
inline Q fma(Q a, Q b, Q c) { // exact when a*b + c cancels (geoDet)
	auto hi = [](Q v) { Q t = v * 0x1p57 + v; return t - (t - v); };
	Q p = a * b, A = hi(a), B = hi(b), a2 = a - A, b2 = b - B;
	return (p + c) + (((A * B - p) + A * b2 + a2 * B) + a2 * b2); }
#endif
}
#define GEO_FN(f) template<class F> F g##f(F x) { \
	if constexpr (is_same_v<F, __float128>) return q128::f(x); \
	else return f(x); }
GEO_FN(sqrt) GEO_FN(sin) GEO_FN(cos) GEO_FN(fabs)
#undef GEO_FN
#define GEO_FN(f) template<class F> F g##f(F x, F y) { \
	if constexpr (is_same_v<F, __float128>) return q128::f(x, y); \
	else return f(x, y); }
GEO_FN(atan2) GEO_FN(hypot)
#undef GEO_FN
template<class F> F gfma(F a, F b, F c) { // a*b+c, rounded once
	if constexpr (is_same_v<F, __float128>) return q128::fma(a, b, c);
	else return fma(a, b, c); }
// a*b - c*d; when the products cancel the rounding errors of both
// are recovered with fma (Kahan), so the sign and digits survive
template<class F> F geoDet(F a, F b, F c, F d) {
	F w = c*d, x = a*b - w;
	if constexpr (geoFloat<F>) {
		if (gfabs(x) > gfabs(w) / 8) return x;
		return gfma(a, b, -w) + gfma(-c, d, w);
	} else return x; }
template<class T> T mulDiv(T a, T b, T c) {
	typedef typename Geo<T>::W W;
	return T(W(a) * W(b) / W(c));
}
template<class T>
struct Point {
	typedef Point P;
	typedef typename Geo<T>::F F;
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
	F dist() const { // hypot only if dist2() over/underflows
		F d = F(dist2()), big = F(sizeof(F) < 8 ? 1e30 : 1e290);
		if (!geoFloat<T> || (d > 1 / big && d < big) ||
			(x == 0 && y == 0)) return gsqrt(d);
		return ghypot(F(x), F(y));
	}
	// angle to x-axis in interval [-pi, pi]
	F angle() const { return gatan2(F(y), F(x)); }
	PF unit() const { // makes dist()=1
		F d = dist(); return PF(F(x)/d, F(y)/d); }
	P perp() const { return P(-y, x); } // rotates +90 degrees
	PF normal() const { return perp().unit(); }
	// returns point rotated 'a' radians ccw around the origin
	PF rotate(F a) const {
		F c = gcos(a), s = gsin(a), u = F(x), v = F(y);
		return PF(u*c - v*s, u*s + v*c); }
	friend ostream& operator<<(ostream& os, P p) {
		return os << "(" << p.x << "," << p.y << ")"; }
};
