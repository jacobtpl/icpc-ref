/**
 * Author: Ulf Lundstrom
 * Date: 2009-02-26
 * License: CC0
 * Source: My head with inspiration from tinyKACTL
 * Description: Class to handle points in the plane, plus the one
 *  block that decides which number types all geometry headers use.
 *  T can be double, long double, __float128, int, long long, __int128
 *  or Fraction. The other geometry headers take their number types
 *  from here (the *-benq.h hull variants have their own `using T`).
 *
 *  Geo<T> (customize below) gives, for coordinate type T:
 *   F: floating type of lengths, angles, ratios and non-integer results.
 *   W: type of intermediates that need more room than T. It must hold
 *      the product of any two values that fit in T.
 *  Defaults: floating T and Fraction compute everything in T itself
 *  (no hidden long double); int -> F=double, W=long long;
 *  long long -> F=double, W=__int128; __int128 -> F=long double, W=__int128
 *  (nothing is wider, so keep values below 2^63 there).
 *
 *  To change precision, change the type you declare your points with:
 *   typedef Point<double> P;       // fast default
 *   typedef Point<long double> P;  // ~19 digits, somewhat slower
 *   typedef Point<__float128> P;   // ~34 digits, much slower
 *  To keep double points but compute circle/half-plane intermediates
 *  in long double (results exact for integer-valued doubles up to 1e9
 *  in the headers that say so), uncomment the Geo<double> line below.
 *  For long double lengths from integer points, edit F in Geo<ll>.
 *
 *  Functions returning non-integer values (intersection points,
 *  centers, areas of unions) return GeoR<T>: F for integer T, else T.
 *
 *  dot, cross and dist2 are evaluated in T: with long long keep
 *  |x|,|y| <= 1e9, with int <= 1.6e4 (prefer long long), with
 *  __int128 <= 4e18. operator/ truncates for integer T.
 *  dist, angle, unit, normal and rotate return F and never truncate.
 *  mulDiv(a,b,c) is a*b/c with the product formed in W (truncates
 *  toward zero for integers).
 *
 *  Epsilon: sgn(x) is the sign test used by the predicates of the
 *  geometry headers. For floating x it treats |x| <= geoEps<type of x>
 *  as zero; the default 0 means exact comparisons. Integers and
 *  Fraction are always exact. Set it once, e.g. in main:
 *   geoEps<double> = 1e-9;
 *  (and geoEps<long double> too if W is long double).
 *  How to choose it: most tested values are cross or dot products of
 *  coordinate differences. With differences up to M their rounding
 *  error is a few u*M*M, where u = 2.2e-16 (double), 1.1e-19
 *  (long double), 1.9e-34 (__float128). eps must be well above that
 *  and below the smallest value that is truly nonzero (for input
 *  with k decimals: 10^(-2k)).
 *   double, M=1e3, 2 decimals: noise 2e-10, smallest 1e-4: use 1e-7.
 *   double, M=1e6, 2 decimals: noise 2e-4 > 1e-4: no eps works;
 *    use long double (noise 1e-7, eps 1e-5) or scale to integers.
 *  If the input is integer-valued, use an integer T and no eps.
 *
 *  __float128 needs no -lquadmath: gsqrt, gsin, gcos, gfabs, gatan2,
 *  ghypot and gfma call glibc's sqrtf128 etc. Read it as long double or
 *  string and print with (long double)x; there is no iostream support.
 * Usage:
 *  typedef Point<double> P;
 *  P a(1, 2), b(4, 6);
 *  double d = (b - a).dist();      // 5
 *  auto c = a.cross(b);            // 1*6 - 2*4 = -2, type T
 *  Point<ll> q(3, 4); double l = q.dist(); // F of long long is double
 * Status: Works fine, used a lot
 */
#pragma once

// ---- Customize below: the only place geometry picks number types ----
template<class T> struct Geo { // floating types and Fraction
	typedef conditional_t<is_class_v<T>, double, T> F; typedef T W; };
template<> struct Geo<int> { typedef double F; typedef ll W; };
template<> struct Geo<ll> { typedef double F; typedef __int128 W; };
template<> struct Geo<__int128> {
	typedef long double F; typedef __int128 W; };
// template<> struct Geo<double> { // wider intermediates for doubles
//	typedef double F; typedef long double W; };
template<class T> T geoEps = T(0); // tolerance of sgn(); 0 = exact
// ---- End of customization ----

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
#define GEO_FN(f) template<class F> F g##f(F x) { \
	if constexpr (is_same_v<F, __float128>) return f##f128(x); \
	else return f(x); }
GEO_FN(sqrt) GEO_FN(sin) GEO_FN(cos) GEO_FN(fabs)
#undef GEO_FN
template<class F> F gatan2(F y, F x) {
	if constexpr (is_same_v<F, __float128>) return atan2f128(y, x);
	else return atan2(y, x); }
template<class F> F ghypot(F x, F y) {
	if constexpr (is_same_v<F, __float128>) return hypotf128(x, y);
	else return hypot(x, y); }
template<class F> F gfma(F a, F b, F c) { // a*b+c, rounded once
	if constexpr (is_same_v<F, __float128>) return fmaf128(a, b, c);
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
