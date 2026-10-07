// Tests 3dHull.h (KACTL) and the Benq hulls: 3dHull-template-benq.h +
// 3dHull-slow-benq.h + 3dHull-fast-benq.h, and the combined 3dHull-benq.h.
// Oracle: exact integer brute force (supporting planes of all triples,
// 2D hull of the points on each plane) giving 6*volume, plus an exact
// certificate check of the returned faces (outward, closed 2-manifold).
#include "../utilities/template.h"
#define mp make_pair
#define pb push_back

#include "../../content/geometry/Point3D.h"
namespace kactl {
#include "../../content/geometry/3dHull.h"
}
#undef E
#undef C
namespace benq {
#include "../../content/geometry/3dHull-template-benq.h"
#include "../../content/geometry/3dHull-slow-benq.h"
#include "../../content/geometry/3dHull-fast-benq.h"
}
namespace benq2 {
#include "../../content/geometry/3dHull-benq.h"
}

typedef __int128 LL;
typedef array<ll,3> I3;
typedef Point3D<double> PD;
typedef array<int,3> Tri;

I3 sub(I3 a, I3 b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
I3 crs(I3 a, I3 b) { return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]}; }
LL dt(I3 a, I3 b) { return (LL)a[0]*b[0] + (LL)a[1]*b[1] + (LL)a[2]*b[2]; }
LL det(I3 a, I3 b, I3 c) { return dt(crs(a, b), c); }
I3 toI(PD p) { return {(ll)p.x, (ll)p.y, (ll)p.z}; }

// 6 * volume of the convex hull, exact; -1 if all points are coplanar
ll bruteVol6(const vector<I3>& p) {
	int n = sz(p);
	set<array<ll,4>> seen;
	LL vol = 0; bool any = 0;
	rep(i,0,n) rep(j,i+1,n) rep(k,j+1,n) {
		I3 nr = crs(sub(p[j], p[i]), sub(p[k], p[i]));
		if (!nr[0] && !nr[1] && !nr[2]) continue;
		bool pos = 0, neg = 0;
		rep(l,0,n) {
			LL d = dt(nr, sub(p[l], p[i]));
			if (d > 0) pos = 1;
			if (d < 0) neg = 1;
		}
		if (pos && neg) continue;
		if (!pos && !neg) return -1;
		if (pos) rep(c,0,3) nr[c] = -nr[c];
		ll g = __gcd(__gcd(abs(nr[0]), abs(nr[1])), abs(nr[2]));
		rep(c,0,3) nr[c] /= g;
		ll off = (ll)dt(nr, p[i]);
		if (!seen.insert({nr[0], nr[1], nr[2], off}).second) continue;
		any = 1;
		int ax = 0;
		rep(c,1,3) if (abs(nr[c]) > abs(nr[ax])) ax = c;
		int u = (ax + 1) % 3, v = (ax + 2) % 3;
		vector<pair<pair<ll,ll>, int>> q;
		rep(l,0,n) if (dt(nr, sub(p[l], p[i])) == 0)
			q.push_back({{p[l][u], p[l][v]}, l});
		sort(all(q));
		q.erase(unique(all(q), [](auto& a, auto& b) { return a.first == b.first; }), q.end());
		auto cr = [&](int a, int b, int c) { // a, b, c index q
			return (LL)(q[b].first.first - q[a].first.first) * (q[c].first.second - q[a].first.second)
			     - (LL)(q[b].first.second - q[a].first.second) * (q[c].first.first - q[a].first.first);
		};
		vi h; // monotone chain, indices into q
		rep(a,0,sz(q)) {
			while (sz(h) >= 2 && cr(h[sz(h)-2], h.back(), a) <= 0) h.pop_back();
			h.push_back(a);
		}
		int lo = sz(h);
		for (int a = sz(q) - 2; a >= 0; a--) {
			while (sz(h) > lo && cr(h[sz(h)-2], h.back(), a) <= 0) h.pop_back();
			h.push_back(a);
		}
		h.pop_back();
		int t = sz(h);
		assert(t >= 3);
		LL part = 0;
		rep(a,1,t-1) part += det(p[q[h[0]].second], p[q[h[a]].second], p[q[h[a+1]].second]);
		// orient outward
		I3 fn = crs(sub(p[q[h[1]].second], p[q[h[0]].second]), sub(p[q[h[2]].second], p[q[h[0]].second]));
		assert(dt(fn, nr) != 0);
		vol += dt(fn, nr) > 0 ? part : -part;
	}
	if (!any) return -1; // all collinear
	return (ll)vol;
}

// exact certificate; returns "" if `fs` is the outward-oriented boundary
// of conv(p) (as a closed triangulated surface with 6*volume = vol6)
string certify(const vector<I3>& p, const vector<Tri>& fs, ll vol6) {
	int n = sz(p);
	set<pii> ed;
	LL vol = 0;
	for (auto& f : fs) {
		rep(c,0,3) if (f[c] < 0 || f[c] >= n) return "index out of range";
		I3 a = p[f[0]], b = p[f[1]], c = p[f[2]];
		I3 nr = crs(sub(b, a), sub(c, a));
		if (!nr[0] && !nr[1] && !nr[2]) return "degenerate face";
		rep(l,0,n) if (dt(nr, sub(p[l], a)) > 0) return "point outside face / face not outward";
		rep(c,0,3) if (!ed.insert({f[c], f[(c+1)%3]}).second) return "duplicate directed edge";
		vol += det(a, b, c);
	}
	for (auto& e : ed) if (!ed.count({e.second, e.first})) return "surface not closed";
	if (vol != vol6) return "wrong volume";
	return "";
}

mt19937_64 gen(12345);
ll rnd(ll lo, ll hi) { return (ll)(gen() % (unsigned long long)(hi - lo + 1)) + lo; }

bool generalPosition(const vector<I3>& p) {
	int n = sz(p);
	rep(i,0,n) rep(j,i+1,n) rep(k,j+1,n) rep(l,k+1,n)
		if (det(sub(p[j], p[i]), sub(p[k], p[i]), sub(p[l], p[i])) == 0) return 0;
	return 1;
}

vector<PD> toD(const vector<I3>& p) {
	vector<PD> r;
	for (auto& q : p) r.push_back(PD((double)q[0], (double)q[1], (double)q[2]));
	return r;
}
vector<I3> toIs(const vector<PD>& p) {
	vector<I3> r;
	for (auto& q : p) r.push_back(toI(q));
	return r;
}

void fail(const char* who, string msg, const vector<I3>& p) {
	cerr << who << ": " << msg << " on n=" << sz(p) << ":";
	for (auto& q : p) cerr << " (" << q[0] << "," << q[1] << "," << q[2] << ")";
	cerr << endl;
	abort();
}

// which: bitmask 1 = kactl, 2 = benq slow, 4 = benq fast
void testAll(const vector<I3>& p, int which, ll vol6 = -2) {
	if (vol6 == -2) vol6 = bruteVol6(p);
	assert(vol6 > 0);
	string r;
	if (which & 1) {
		vector<PD> d = toD(p);
		vector<Tri> fs;
		for (auto& f : kactl::hull3d(d)) fs.push_back({f.a, f.b, f.c});
		if ((r = certify(p, fs, vol6)) != "") fail("3dHull.h", r, p);
	}
	if (which & 2) {
		// the Benq hulls permute their argument; faces index the permuted vector
		vector<PD> d = toD(p);
		vector<Tri> fs = benq::hull3d(d);
		if ((r = certify(toIs(d), fs, vol6)) != "") fail("3dHull-slow-benq.h", r, p);
		d = toD(p); fs = benq2::hull3d(d);
		if ((r = certify(toIs(d), fs, vol6)) != "") fail("3dHull-benq.h hull3d", r, p);
	}
	if (which & 4) {
		vector<PD> d = toD(p);
		vector<Tri> fs = benq::hull3dFast(d);
		if ((r = certify(toIs(d), fs, vol6)) != "") fail("3dHull-fast-benq.h", r, p);
		d = toD(p); fs = benq2::hull3dFast(d);
		if ((r = certify(toIs(d), fs, vol6)) != "") fail("3dHull-benq.h hull3dFast", r, p);
	}
}

// larger inputs: no brute force, certify with the volume of one hull and
// require every implementation to agree
ll hullVol6(const vector<I3>& p, const vector<Tri>& fs) {
	LL v = 0;
	for (auto& f : fs) v += det(p[f[0]], p[f[1]], p[f[2]]);
	return (ll)v;
}

int main() {
	// 1. general position (the documented precondition), tiny to small n
	for (int C : {3, 10, 1000, 10000}) rep(it,0,6000) {
		int n = (int)rnd(4, C == 3 ? 7 : 10);
		vector<I3> p(n);
		do {
			for (auto& q : p) rep(c,0,3) q[c] = rnd(-C, C);
		} while (!generalPosition(p));
		testAll(p, 7);
	}
	// 2. general position, points on a sphere-ish shell (all on the hull),
	// sorted input
	rep(it,0,300) {
		int n = (int)rnd(4, 14);
		vector<I3> p(n);
		do {
			for (auto& q : p) {
				double a = (double)rnd(0, 1000000) / 1e6 * 6.283185307, z = (double)rnd(-1000000, 1000000) / 1e6;
				double r = sqrt(1 - z * z);
				q = {(ll)(1000 * r * cos(a)), (ll)(1000 * r * sin(a)), (ll)(1000 * z)};
			}
		} while (!generalPosition(p));
		sort(all(p));
		testAll(p, 7);
	}
	// 3. degenerate inputs: coplanar / collinear / duplicate points (but not
	// all coplanar). This violates the precondition documented for 3dHull.h
	// (which does fail here; define HULL3D_KACTL_DEGENERATE to see), but the
	// Benq hulls handle it.
	rep(it,0,30000) {
		int n = (int)rnd(4, it % 10 ? 9 : 40), C = (int)rnd(1, 3);
		vector<I3> p(n);
		ll v;
		do {
			for (auto& q : p) rep(c,0,3) q[c] = rnd(0, C);
		} while ((v = bruteVol6(p)) <= 0);
#ifdef HULL3D_KACTL_DEGENERATE
		testAll(p, 7, v);
#else
		testAll(p, 6, v);
#endif
	}
	// 3b. same, more points (dense small grid: many duplicates and many
	// coplanar points)
	rep(it,0,40) {
		int n = (int)rnd(50, 110), C = (int)rnd(2, 5);
		vector<I3> p(n);
		ll v;
		do {
			for (auto& q : p) rep(c,0,3) q[c] = rnd(0, C);
		} while ((v = bruteVol6(p)) <= 0);
		testAll(p, 6, v);
	}
	// 4. medium n, random points in a cube with large coordinates (four
	// coplanar points are astronomically unlikely); all hulls must agree
	rep(it,0,30) {
		int n = (int)rnd(50, 300);
		vector<I3> p(n);
		for (auto& q : p) rep(c,0,3) q[c] = rnd(-20000, 20000);
		vector<PD> d = toD(p);
		vector<Tri> fs = benq::hull3dFast(d);
		ll v = hullVol6(toIs(d), fs);
		testAll(p, 7, v);
	}
#ifdef HULL3D_EXPECT_NO_REORDER
	// The Benq hulls shuffle their argument and the returned faces index
	// the shuffled vector, not the caller's original order.
	{
		vector<I3> p = {{0,0,0}, {5,1,0}, {1,6,2}, {2,1,7}, {6,7,8}, {-3,2,1}, {1,-4,3}, {2,3,-5}};
		vector<PD> d = toD(p);
		benq::hull3dFast(d);
		assert(toIs(d) == p);
	}
#endif
#ifdef HULL3D_BIG_COPLANAR
	// Coplanar points whose orientation determinants (~coord^3) are not
	// exact in double: the Benq hulls return an invalid hull.
	testAll({{907442,883669,-971140}, {210506,2544043,-2133836}, {1278316,-558966,371986},
		{1837290,1154890,-1180502}, {907442,883669,-971140}, {1837290,1154890,-1180502}}, 6);
#endif
	cout << "Tests passed!" << endl;
}
