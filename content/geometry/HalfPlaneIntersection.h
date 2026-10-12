/**
 * Author: Benq
 * Description: Returns vertices of half-plane intersection in CCW order.
   * A half-plane is the area to the left of a ray, which is defined
   * by a point \texttt{p} and a nonzero direction \texttt{dp}.
   * Returns \{\} unless the intersection is bounded with positive area
   * (so also when it is empty, a point, a segment or unbounded).
   * With \texttt{RayT<ll>} every decision is exact (\_\_int128 products) for
   * $|p| \le 10^9$, $|dp| \le 2 \cdot 10^9$ and only the returned vertices
   * are rounded (long double); use it for integer input.
   * Floating rays with integer-valued coordinates within these bounds
   * also use the exact integer predicates.
   * With floating point rays products are formed in long double and a
   * determinant within \texttt{eps} times the sum of the magnitudes of
   * its terms counts as zero, i.e. the configuration as degenerate
   * (parallel rays, or a vertex on a third line). Keep \texttt{eps} near
   * rounding error: features smaller than \texttt{eps} times the distance
   * between the origins \texttt{p} are lost. The tests are invariant under
   * translation and under scaling of the input or of any \texttt{dp}.
   * \texttt{add\_bounds} clips to $[0,10^9]^2$ (not to negative coordinates).
 * Time: O(N\log N)
 * Source: Own
 * Status: stress-tested against exact brute force; older version tested on
	* https://open.kattis.com/problems/bigbrother 
	* https://icpc.kattis.com/problems/domes
	* https://codeforces.com/gym/100962 B
 */
#pragma once

#include "Point.h"

template<class T> struct RayT { Point<T> p, dp; };
typedef RayT<double> Ray;

template<class T, class F =
	conditional_t<is_integral_v<T>, long double, T>>
vector<Point<F>> halfPlaneIsect(vector<RayT<T>> rays,
	bool add_bounds = false, long double eps = 1e-13) {
	typedef Point<T> P; typedef RayT<T> R;
	if constexpr (is_floating_point_v<T>) {
		auto integer = [](T x, T bound) {
			return abs(x) <= bound && trunc(x) == x; };
		if (all_of(all(rays), [&](const R& r) {
			return integer(r.p.x, T(1e9)) && integer(r.p.y, T(1e9)) &&
				integer(r.dp.x, T(2e9)) && integer(r.dp.y, T(2e9)); })) {
			vector<RayT<ll>> exact;
			for (const R& r : rays) exact.push_back({
				Point<ll>((ll)r.p.x, (ll)r.p.y), Point<ll>((ll)r.dp.x, (ll)r.dp.y)});
			vector<Point<F>> res;
			for (auto p : halfPlaneIsect(exact, add_bounds))
				res.emplace_back(F(p.x), F(p.y));
			return res;
		}
	}
	const bool I = is_integral_v<T>;
	typedef conditional_t<is_integral_v<T>, __int128, long double> W;
	W E = I ? W(0) : W(eps);
	if (add_bounds) { // bound input by rectangle [0,D] x [0,D]
		T D = (T)1e9;
		rays.push_back({P(0,0),P(1,0)});
		rays.push_back({P(D,0),P(0,1)});
		rays.push_back({P(D,D),P(-1,0)});
		rays.push_back({P(0,D),P(0,-1)});
	}
	auto cr = [](P a, P b) { return (W)a.x*(W)b.y - (W)a.y*(W)b.x; };
	auto ab = [](P a, P b) { // cr with absolute values of the terms
		W x = (W)a.x * (W)b.y, y = (W)a.y * (W)b.x;
		return (x < 0 ? -x : x) + (y < 0 ? -y : y); };
	// -1 if lower half or +x axis, 1 if upper half or -x axis
	auto half = [](P x) { return x.y != 0 ? sgn(x.y) : -sgn(x.x); };
	sort(all(rays), [&](const R& a, const R& b) {
		int A = half(a.dp), B = half(b.dp);
		return A == B ? cr(a.dp, b.dp) > 0 : A < B; });
	// 1: b is left of a by less than half a turn, 0: same direction
	auto turn = [&](P a, P b) {
		W x = cr(a, b), y = (W)a.x * (W)b.x + (W)a.y * (W)b.y;
		W t = E * (y < 0 ? -y : y);
		return x > t ? 1 : x < -t || y <= 0 ? -1 : 0;
	};
	// is isect(a,b) not strictly left of c? needs turn(a.dp,b.dp) = 1
	auto out = [&](const R& a, const R& b, const R& c) {
		P ca = a.p - c.p, ab_ = b.p - a.p;
		W u = cr(c.dp, ca) * cr(a.dp, b.dp)
			+ cr(ab_, b.dp) * cr(c.dp, a.dp);
		return u <= (I ? W(0) : E * (ab(c.dp, ca) * ab(a.dp, b.dp)
			+ ab(ab_, b.dp) * ab(c.dp, a.dp)));
	};
	vector<R> d; // rays are only removed if redundant; returning {}
	for (auto& c : rays) { // means no interior or a turn >= 180
		if (sz(d) && !turn(d.back().dp, c.dp)) { // keep the left one
			if (cr(c.dp, c.p - d.back().p) <= 0) continue;
			d.pop_back();
		}
		while (sz(d)) {
			if (turn(d.back().dp, c.dp) < 1) return {};
			if (sz(d) < 2 || !out(d[sz(d)-2], d.back(), c)) break;
			d.pop_back();
		}
		d.push_back(c);
	}
	int h = 0, t = sz(d); // answer is d[h..t), cyclically
	while (t - h > 2) {
		if (turn(d[t-1].dp, d[h].dp) < 1) return {};
		if (out(d[t-2], d[t-1], d[h])) t--;
		else if (out(d[h], d[h+1], d[t-1])) h++;
		else break;
	}
	vector<Point<F>> res;
	if (t - h > 2) rep(i,h,t) {
		const R &a = d[i], &b = d[i + 1 == t ? h : i + 1];
		F s = F(cr(b.p - a.p, b.dp)) / F(cr(a.dp, b.dp));
		res.push_back(Point<F>(F(a.p.x) + F(a.dp.x) * s,
			F(a.p.y) + F(a.dp.y) * s));
	}
	return res;
}
