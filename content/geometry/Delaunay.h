/**
 * Author: Benq, adapted
 * Date: Unknown
 * Description: Computes the Delaunay triangulation of a set of points
 * as a list of counter-clockwise triangles.
 * Each circumcircle contains none of the input points.
 * Wrapper around FastDelaunay.
 * Duplicate points are merged. Fewer than three distinct points or
 * all points on a line give no triangles; cocircular points give
 * any valid triangulation.
 * Exact for integer coordinates of absolute value up to $10^9$.
 * Time: O(n \log n)
 * Status: stress-tested against exact oracle
 */
#pragma once

#include "FastDelaunay.h"

vector<array<P,3>> triHull(vector<P> p) {
	sort(all(p)); p.erase(unique(all(p)), p.end());
	vector<P> t = triangulate(p);
	vector<array<P,3>> res;
	for (int i = 0; i < sz(t); i += 3)
		res.push_back({t[i], t[i+1], t[i+2]});
	return res;
}
