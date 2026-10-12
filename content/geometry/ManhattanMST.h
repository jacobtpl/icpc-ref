/**
 * Author: chilli, Takanori MAEHARA
 * Date: 2019-11-02
 * License: CC0
 * Source: https://github.com/spaghetti-source/algorithm/blob/master/geometry/rectilinear_mst.cc
 * Description: Given N points, returns up to 4*N edges, which are guaranteed
 * to contain a minimum spanning tree for the graph with edge weights w(p, q) =
 * |p.x - q.x| + |p.y - q.y|. Edges are in the form (distance, src, dst). Use a
 * standard MST algorithm on the result to find the final MST.
 * Takes Point<int> or Point<ll>; all arithmetic and the returned weights are
 * ll. Exact for any int input, and for ll input with $|x|, |y| \le 2 \cdot 10^{18}$
 * (weights are at most 4 max|coord|). The MST total is the caller's sum:
 * it fits in ll while $(N-1) \cdot 4 \max|coord| < 2^{63}$.
 * Time: O(N \log N)
 * Status: Stress-tested
 */
#pragma once
#include "Point.h"

template<class T>
vector<array<ll, 3>> manhattanMST(const vector<Point<T>>& in) {
	int n = sz(in);
	vector<array<ll, 2>> ps(n);
	rep(i,0,n) ps[i] = {in[i].x, in[i].y};
	vi id(n);
	iota(all(id), 0);
	vector<array<ll, 3>> edges;
	rep(k,0,4) {
		sort(all(id), [&](int i, int j) {
		     return ps[i][0] + ps[i][1] < ps[j][0] + ps[j][1];});
		map<ll, int> sweep;
		for (int i : id) {
			for (auto it = sweep.lower_bound(-ps[i][1]);
				        it != sweep.end(); sweep.erase(it++)) {
				int j = it->second;
				ll dx = ps[i][0] - ps[j][0], dy = ps[i][1] - ps[j][1];
				if (dy > dx) break;
				edges.push_back({dx + dy, i, j});
			}
			sweep[-ps[i][1]] = i;
		}
		for (auto& p : ps) if (k & 1) p[0] = -p[0]; else swap(p[0], p[1]);
	}
	return edges;
}
