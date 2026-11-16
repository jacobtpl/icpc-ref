/**
 * Author: chilli, Takanori MAEHARA
 * Date: 2019-11-02
 * License: CC0
 * Source: https://github.com/spaghetti-source/algorithm/blob/master/geometry/rectilinear_mst.cc
 * Description: Given N points, returns up to 4*N edges, which are guaranteed
 * to contain a minimum spanning tree for the graph with edge weights w(p, q) =
 * |p.x - q.x| + |p.y - q.y|. Edges are in the form (distance, src, dst). Use a
 * standard MST algorithm on the result to find the final MST.
 * Takes any Point<T>. All arithmetic and the returned weights use D,
 * which defaults to T but at least long long (int -> ll, ll -> ll,
 * __int128 -> __int128, double -> double, ...); override it with
 * manhattanMST<T, D>(pts). Vertex indices are stored as D too.
 * Weights are at most 4 max|coord|: exact for any int input, and for
 * ll input with |x|, |y| <= 2e18. The MST total is the caller's sum.
 * Usage:
 *  vector<Point<int>> pts = ...;
 *  for (auto [w, i, j] : manhattanMST(pts)) ...; // w, i, j are ll
 * Time: O(N log N)
 * Status: Stress-tested
 */
#pragma once
#include "Point.h"

template<class T, class D = common_type_t<T, ll>>
vector<array<D, 3>> manhattanMST(const vector<Point<T>>& in) {
	int n = sz(in);
	vector<array<D, 2>> ps(n);
	rep(i,0,n) ps[i] = {D(in[i].x), D(in[i].y)};
	vi id(n);
	iota(all(id), 0);
	vector<array<D, 3>> edges;
	rep(k,0,4) {
		sort(all(id), [&](int i, int j) {
		     return ps[i][0] + ps[i][1] < ps[j][0] + ps[j][1];});
		map<D, int> sweep;
		for (int i : id) {
			for (auto it = sweep.lower_bound(-ps[i][1]);
				        it != sweep.end(); sweep.erase(it++)) {
				int j = it->second;
				D dx = ps[i][0] - ps[j][0], dy = ps[i][1] - ps[j][1];
				if (dy > dx) break;
				edges.push_back({dx + dy, D(i), D(j)});
			}
			sweep[-ps[i][1]] = i;
		}
		for (auto& p : ps) if (k & 1) p[0] = -p[0]; else swap(p[0], p[1]);
	}
	return edges;
}
