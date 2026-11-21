/**
 * Author: Lukas Polacek
 * Date: 2009-10-28
 * License: CC0
 * Source:
 * Description: Simple bipartite matching (Kuhn). Same interface as
 *  hopcroftKarp.h: g[a] lists the right neighbours of left vertex a, btoa has
 *  one entry per right vertex, all -1. Returns the matching size; btoa[b]
 *  becomes the left vertex matched to b, or -1.
 *  Left vertices are tried in order 0..n-1 and a matched left vertex stays
 *  matched, so sorting the left side gives a greedy / weighted-by-left match.
 * Time: O(VE)
 * Usage: vi btoa(m, -1); dfsMatching(g, btoa);
 * Status: stress-tested by MinimumVertexCover
 */
#pragma once

int dfsMatching(vector<vi>& g, vi& btoa) {
	vi vis;
	auto go = [&](auto& self, int j) -> bool {
		if (btoa[j] == -1) return 1;
		vis[j] = 1; int di = btoa[j];
		for (int e : g[di])
			if (!vis[e] && self(self, e)) {
				btoa[e] = di;
				return 1;
			}
		return 0;
	};
	rep(i,0,sz(g)) {
		vis.assign(sz(btoa), 0);
		for (int j : g[i])
			if (go(go, j)) {
				btoa[j] = i;
				break;
			}
	}
	return sz(btoa) - (int)count(all(btoa), -1);
}
