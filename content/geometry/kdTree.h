/**
 * Author: Stanford
 * Date: Unknown
 * Source: Stanford Notebook
 * Description: KD-tree (2d, can be extended to 3d). Nodes live in a vector
 * owned by the tree and refer to each other by index, so trees are freed on
 * destruction and can be copied, moved and assigned freely.
 * KDTree<P> works on any P = Point<T>. Squared distances are
 * evaluated in T: with long long, tree points and queries need
 * |x|, |y| <= 2^30-1 (about 1.07e9), so that
 * dist2 <= 8 * (2^30-1)^2 < 2^63. Duplicate points are fine.
 * nearest on an empty tree returns {-1, P()}.
 * Usage:
 *  typedef Point<ll> P;
 *  KDTree<P> tree(pts);
 *  auto [d2, q] = tree.nearest(P(3, 4)); // squared distance, point
 * Time: build O(N log N). nearest is O(log N) on average for random
 * points, but O(N) in the worst case (e.g. points on a circle or diagonal).
 * Status: Stress-tested
 */
#pragma once

#include "Point.h"

template<class P> struct KDTree {
	typedef decltype(P().x) T;
	static bool on_x(const P& a, const P& b) { return a.x < b.x; }
	static bool on_y(const P& a, const P& b) { return a.y < b.y; }

	struct Node {
		P pt; // if this is a leaf, the single point in it
		T x0, x1, y0, y1; // bounds
		int first = -1, second = -1; // children (indices), -1 if leaf

		T distance(const P& p) const { // min squared distance to a point
			T x = (p.x < x0 ? x0 : p.x > x1 ? x1 : p.x);
			T y = (p.y < y0 ? y0 : p.y > y1 ? y1 : p.y);
			return (P(x,y) - p).dist2();
		}
	};

	vector<Node> t; // t[0] is the root (if any)

	int build(vector<P>& vp, int l, int r) { // builds from vp[l, r)
		int id = sz(t);
		t.emplace_back();
		Node n;
		n.pt = vp[l];
		n.x0 = n.x1 = vp[l].x; n.y0 = n.y1 = vp[l].y;
		for (int i = l; i < r; i++) {
			n.x0 = min(n.x0, vp[i].x); n.x1 = max(n.x1, vp[i].x);
			n.y0 = min(n.y0, vp[i].y); n.y1 = max(n.y1, vp[i].y);
		}
		if (r - l > 1) {
			// split on x if width >= height (not ideal...)
			// divide by taking half the array for each child (not
			// best performance with many duplicates in the middle)
			int half = l + (r - l) / 2;
			nth_element(vp.begin() + l, vp.begin() + half,
				vp.begin() + r, n.x1 - n.x0 >= n.y1 - n.y0 ? on_x : on_y);
			n.first = build(vp, l, half);
			n.second = build(vp, half, r);
		}
		t[id] = n;
		return id;
	}

	KDTree(vector<P> vp = {}) {
		if (vp.empty()) return;
		t.reserve(2 * vp.size() - 1);
		build(vp, 0, sz(vp));
	}

	pair<T, P> search(int i, const P& p) const {
		const Node& node = t[i];
		if (node.first < 0) {
			// uncomment if we should not find the point itself:
			// if (p == node.pt) return {infinity, P()};
			return make_pair((p - node.pt).dist2(), node.pt);
		}

		int f = node.first, s = node.second;
		T bfirst = t[f].distance(p), bsec = t[s].distance(p);
		if (bfirst > bsec) swap(bsec, bfirst), swap(f, s);

		// search closest side first, other side if needed
		auto best = search(f, p);
		if (bsec < best.first)
			best = min(best, search(s, p));
		return best;
	}

	// find nearest point to a point, and its squared distance
	// (requires an arbitrary operator< for Point)
	pair<T, P> nearest(const P& p) const {
		return t.empty() ? make_pair(T(-1), P()) : search(0, p);
	}
};
