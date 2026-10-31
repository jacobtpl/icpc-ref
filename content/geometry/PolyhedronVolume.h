/**
 * Author: Mattias de Zalenski
 * Date: 2002-11-04
 * Description: Magic formula for the volume of a polyhedron. Faces should point outwards.
 * Status: tested
 */
#pragma once

template<class V, class L>
auto signedPolyVolume(const V& p, const L& trilist) {
	using T = decltype(p[0].dot(p[0]));
	conditional_t<is_integral<T>::value, double, T> v = 0;
	for (auto i : trilist) v += p[i.a].cross(p[i.b]).dot(p[i.c]);
	return v / 6;
}
