/**
 * Author: Mattias de Zalenski
 * Date: 2002-11-04
 * Description: Magic formula for the volume of a polyhedron. Faces should point outwards.
 *  Faces are triangles given as \{a,b,c\} members or as array<int,3> (both hull versions).
 *  Coordinates are taken relative to a vertex, so the result does not degrade
 *  far from the origin. signedPolyVolume6 returns six times the volume in the
 *  coordinate type: exact for integers while (number of faces) $\cdot 6D^3$ fits,
 *  $D$ being the largest coordinate difference. signedPolyVolume returns a
 *  double for integer coordinates (also \_\_int128) and otherwise the coordinate type (Fraction is exact).
 * Status: tested
 */
#pragma once

template<class F> auto triIdx(const F& f, int) ->
	decltype(array<int,3>{{f.a, f.b, f.c}}) { return {{f.a, f.b, f.c}}; }
template<class F> auto triIdx(const F& f, long) ->
	decltype(array<int,3>{{f[0], f[1], f[2]}}) {
	return {{f[0], f[1], f[2]}}; }

template<class V, class L>
auto signedPolyVolume6(const V& p, const L& trilist) {
	decltype(p[0].x) v = 0;
	int o = -1;
	for (auto& f : trilist) {
		auto i = triIdx(f, 0);
		if (o < 0) o = i[0];
		v += (p[i[0]]-p[o]).cross(p[i[1]]-p[o]).dot(p[i[2]]-p[o]);
	}
	return v;
}
template<class V, class L>
auto signedPolyVolume(const V& p, const L& trilist) {
	typedef decltype(p[0].x) T;
	typedef conditional_t<is_floating_point<T>::value ||
		is_class<T>::value, T, double> R; // ll, __int128 -> double
	return R(signedPolyVolume6(p, trilist)) / 6;
}
