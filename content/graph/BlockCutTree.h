/**
 * Author: Siyong
 * Date: 2024-04-07
 * License: CC0
 * Source: folklore
 * Description: Builds the block cut tree.
 * 	BCTree node $n$ is an AP if $n >= cut$. Corresponds to who[n][0] in original graph
 * 	Node $v$ is an AP if $vmap[v] >= cut$.
 * 	$emap[i] = -1$ if edge $i$ is a bridge. Otherwise, $emap[i] = $ the BCC containing it
 * Usage:
 *  see BiconnectedComponents.h
 *  edges[i] = the edge i (pair of two nodes)
 * Time: O(E + V)
 * Status: tested during MIPT ICPC Workshop 2017
 */
#pragma once

#include "BiconnectedComponents.h"

vector<pii> edges;
tuple<int, vector<vi>, vector<vi>, vi, vi> BCTree() {
	int N = ed.size(), M = edges.size();
	vector<int> emap(M, -1); // edge -> bicomp id
	vector<vi> bclist; // list of biconnected components
	bicomps([&](vector<int> &&eds) {
		for(int x: eds) emap[x] = bclist.size();
		bclist.emplace_back(eds);
	});

	vector<int> vmap(N, -1);
	for(int i = 0;i < M; ++i)
		if(emap[i] == -1) { // bridge: connects two APs
			auto [u, v] = edges[i];
			vmap[u] = vmap[v] = -2;
		}

	for(int i = 0;i < bclist.size();++i)
		for(int x: bclist[i]) {
			auto [u, v] = edges[x];
			for(int j = 2;j--; swap(u, v))
				if(vmap[u] == -1) vmap[u] = i;
				else if(vmap[u] != i) vmap[u] = -2;
		}

	int const cut = bclist.size();
	int TN = bclist.size();
	vector<vi> who(TN);
	for(int i = 0;i < N;++i)
		if(vmap[i] == -2) vmap[i] = TN++, who.emplace_back(1, i);
		else who[vmap[i]].emplace_back(i);
	vector<vi> tadj(TN);

	for(int i = 0;i < N; ++i)
		if(cut <= vmap[i]) // if `i` is an AP
			for(auto [x, e]: ed[i]) {
				if(emap[e] == -1) // Bridge: connect both APs
					tadj[vmap[i]].push_back(vmap[x]);
				else {
					tadj[vmap[i]].push_back(emap[e]);
					tadj[emap[e]].push_back(vmap[i]);
				}
			}
	for(auto &v: tadj) { // one AP can connect to a BCC in multiple ways
		sort(all(v));
		v.resize(distance(v.begin(), unique(all(v))));
	}
	return {cut, tadj, who, emap, vmap};
}
