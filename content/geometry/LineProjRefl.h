/**
 * Author: Victor Lecomte, chilli
 * Date: 2019-10-29
 * License: CC0
 * Source: https://vlecomte.github.io/cp-geo.pdf
 * Description: Projects point p onto line ab. Set refl=true to get reflection
 * of point p across line ab instead. Requires a != b (asserted). The result
 * is p plus an offset, so doubles lose no digits to the absolute position.
 * For Point<ll>
 * coordinates must be at most $10^9$ in absolute value (triple products are
 * done in \_\_int128) and each coordinate of the offset is truncated toward
 * zero when the desired point doesn't have integer coordinates.
 * Status: stress-tested
 */
#pragma once

#include "Point.h"

template<class P>
P lineProj(P a, P b, P p, bool refl=false) {
	P v = b - a, w = v.perp() * (1 + refl);
	auto c = v.cross(p - a), d = v.dist2();
	assert(d != 0);
	return p - P(mulDiv(w.x, c, d), mulDiv(w.y, c, d));
}
