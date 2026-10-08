/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: folklore (CDQ divide and conquer, after Chen Danqi); written from scratch
 * Description: CDQ divide and conquer: sort by the first key, recurse on
 *  both halves, then add the contribution of the left half to the right
 *  half with a sweep over the second key and a Fenwick tree over the third.
 *  Shown on 3D partial order: res[i] is the number of $j \neq i$ with
 *  $p_j \le p_i$ in all three coordinates (equal points count each other).
 *  Any int coordinates work (z is compressed). To drop a dimension by
 *  offline processing in general, make time the first key and change the
 *  marked lines.
 * Time: O(N \log^2 N)
 * Usage: vi res = cdq(pts); // pts[i] = \{x, y, z\}
 * Status: stress-tested
 */
#pragma once

#include "../data-structures/FenwickTree.h"

vi cdq(vector<array<int, 3>> p) {
	int n = sz(p);
	vi zs(n), id(n), res(n);
	rep(i,0,n) zs[i] = p[i][2];
	sort(all(zs));
	for (auto& a : p) a[2] = int(lower_bound(all(zs), a[2]) - zs.begin());
	iota(all(id), 0);
	sort(all(id), [&](int a, int b) { return p[a] < p[b]; });
	vi o = id;
	FT ft(n);
	auto cmp = [&](int a, int b) { return p[a][1] < p[b][1]; };
	auto rec = [&](auto& self, int l, int r) -> void {
		if (r - l < 2) return;
		int m = (l + r) / 2, i = l;
		self(self, l, m); self(self, m, r);
		// id[l,m), id[m,r) are sorted by y; all of left <= right in x
		rep(j,m,r) {
			for (; i < m && !cmp(id[j], id[i]); i++)
				ft.update(p[id[i]][2], 1); // insert left
			res[id[j]] += (int)ft.query(p[id[j]][2] + 1); // query right
		}
		rep(k,l,i) ft.update(p[id[k]][2], -1); // undo
		inplace_merge(id.begin()+l, id.begin()+m, id.begin()+r, cmp);
	};
	rec(rec, 0, n);
	for (int i = n - 1; i-- > 0;) // equal points share the last one's answer
		if (p[o[i]] == p[o[i+1]]) res[o[i]] = res[o[i+1]];
	return res;
}
