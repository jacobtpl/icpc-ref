/**
 * Author: Johan Sannemo
 * Date: 2017-04-18
 * Source: derived from https://gist.github.com/msg555/4963794 by Mark Gordon
 * Description: Computes all faces of the 3-dimension hull of a point set.
 *  *No four points must be coplanar*, or else random results will be returned.
 *  The first four points must not be coplanar either.
 *  All faces will point outwards.
 *  Works on any Point3D<T>; all tests use coordinate differences, so
 *  only the extent $D$ (largest coordinate difference) matters:
 *  integer input is exact while $6D^3$ fits in T (double: $D \le 1.1 \cdot 10^5$
 *  if the coordinates themselves are below $2^{53}$, long long: $D \le 10^6$,
 *  \_\_int128: $D \le 2 \cdot 10^{12}$).
 * Time: O(n^2)
 * Status: tested on SPOJ CH3D
 */
#pragma once

#include "Point3D.h"

template<class P3> struct Face { P3 q; int a, b, c; };
typedef Point3D<double> P3;
typedef Face<P3> F;
struct PR {
	void ins(int x) { (a == -1 ? a : b) = x; }
	void rem(int x) { (a == x ? a : b) = -1; }
	int cnt() { return (a != -1) + (b != -1); }
	int a, b;
};

template<class P3>
vector<Face<P3>> hull3d(const vector<P3>& A) {
	typedef Face<P3> F;
	assert(sz(A) >= 4);
	vector<vector<PR>> E(sz(A), vector<PR>(sz(A), {-1, -1}));
#define E(x,y) E[f.x][f.y]
	vector<F> FS;
	auto mf = [&](int i, int j, int k, int l) {
		P3 q = (A[j] - A[i]).cross((A[k] - A[i]));
		if (q.dot(A[l] - A[i]) > 0)
			q = q * -1;
		F f{q, i, j, k};
		E(a,b).ins(k); E(a,c).ins(j); E(b,c).ins(i);
		FS.push_back(f);
	};
	rep(i,0,4) rep(j,i+1,4) rep(k,j+1,4)
		mf(i, j, k, 6 - i - j - k);

	rep(i,4,sz(A)) {
		rep(j,0,sz(FS)) {
			F f = FS[j];
			if(f.q.dot(A[i] - A[f.a]) > 0) {
				E(a,b).rem(f.c);
				E(a,c).rem(f.b);
				E(b,c).rem(f.a);
				swap(FS[j--], FS.back());
				FS.pop_back();
			}
		}
		int nw = sz(FS);
		rep(j,0,nw) {
			F f = FS[j];
#define C(a, b, c) if (E(a,b).cnt() != 2) mf(f.a, f.b, i, f.c);
			C(a, b, c); C(a, c, b); C(b, c, a);
		}
	}
	// orient outwards: no dot of two normals (overflows integers)
	for (F& it : FS) rep(l,0,4) if ((A[it.b] - A[it.a]).cross(
		A[it.c] - A[it.a]).dot(A[l] - A[it.a]) > 0) {
		swap(it.c, it.b); break; }
	return FS;
}
#undef E
#undef C
