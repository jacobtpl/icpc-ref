/**
 * Author: chilli, Takanori MAEHARA
 * Date: 2020-04-03
 * License: CC0
 * Source: https://github.com/spaghetti-source/algorithm/blob/master/graph/gomory_hu_tree.cc#L102
 * Description: Given the edges {a, b, cap} of an undirected flow graph,
 *  returns the edges {a, b, cut} of the Gomory-Hu tree. The min cut between
 *  two vertices is the minimum edge weight on their tree path.
 *  Any flow struct with addEdge(a, b, c, rcap), calc and leftOfMinCut works
 *  (PushRelabel.h has the same interface).
 * Time: V - 1 max flows
 * Status: Tested on CERC 2015 J, stress-tested
 */
#pragma once

#include "Dinic.h"

typedef array<ll, 3> Edge;
vector<Edge> gomoryHu(int N, vector<Edge> ed) {
	vector<Edge> tree;
	vi par(N);
	rep(i,1,N) {
		Dinic D(N);
		for (Edge t : ed) D.addEdge((int)t[0], (int)t[1], t[2], t[2]);
		tree.push_back({i, par[i], D.calc(i, par[i])});
		rep(j,i+1,N)
			if (par[j] == par[i] && D.leftOfMinCut(j)) par[j] = i;
	}
	return tree;
}
