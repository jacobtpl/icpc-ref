// Exact validation of the Benq hulls (3dHull-template/slow/fast-benq.h) in
// their documented domain (integer coordinates, no four coplanar), far from
// the origin, and of PolyhedronVolume.h on their array<int,3> faces.
#include "../utilities/template.h"
#define mp make_pair
#define pb push_back
using pi = pair<int,int>;

#include "../../content/geometry/Point3D.h"
#include "../../content/geometry/3dHull-template-benq.h"
#include "../../content/geometry/3dHull-slow-benq.h"
#include "../../content/geometry/3dHull-fast-benq.h"
#include "../../content/geometry/PolyhedronVolume.h"

typedef __int128 LL;
typedef array<LL,3> IP;

mt19937_64 gen64(20261012);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(gen64); }

LL det(const IP& a, const IP& b, const IP& c, const IP& d) {
	LL u[3], v[3], w[3];
	rep(i,0,3) u[i] = b[i]-a[i], v[i] = c[i]-a[i], w[i] = d[i]-a[i];
	return (u[1]*v[2] - u[2]*v[1]) * w[0] + (u[2]*v[0] - u[0]*v[2]) * w[1]
		+ (u[0]*v[1] - u[1]*v[0]) * w[2];
}
F canon(F t) {
	while (t[0] > t[1] || t[0] > t[2]) rotate(t.begin(), t.begin()+1, t.end());
	return t;
}
void fail(const string& why, const vector<P3>& p) {
	cout << why << "\nPoints:";
	for (auto& q : p) cout << " (" << (double)q.x << "," << (double)q.y
		<< "," << (double)q.z << ")";
	cout << endl;
	abort();
}
// p is in the order left by the hull call (prep() permutes it)
void check(const string& name, const vector<P3>& p, const vector<F>& fs) {
	int n = sz(p);
	vector<IP> q;
	for (auto& a : p) q.push_back({(LL)a.x, (LL)a.y, (LL)a.z});
	set<F> want, got;
	rep(i,0,n) rep(j,i+1,n) rep(k,j+1,n) {
		int pos = 0, neg = 0;
		rep(l,0,n) if (l != i && l != j && l != k)
			(det(q[i], q[j], q[k], q[l]) > 0 ? pos : neg)++;
		if (!pos) want.insert(canon({i, j, k}));
		else if (!neg) want.insert(canon({i, k, j}));
	}
	for (auto& f : fs) if (!got.insert(canon(f)).second)
		fail(name + ": repeated face", p);
	if (got != want) fail(name + ": wrong faces", p);
	LL v = 0;
	for (auto& t : want) v += det(q[0], q[t[0]], q[t[1]], q[t[2]]);
	if ((LL)signedPolyVolume6(p, fs) != v) fail(name + ": volume", p);
	if (abs(signedPolyVolume(p, fs) * 6 - (double)v) > 1e-11 * (double)v)
		fail(name + ": volume (double)", p);
}

int main() {
	int cases = 0;
	// T = double: exact for |coordinate differences| <= 1e5
	rep(it,0,1500) for (ll R : {4LL, 50LL, 1000LL}) rep(slab,0,2)
		for (ll off : {0LL, 1000000000LL, 10000000000000LL,
				4000000000000000LL}) {
			int n = (int)rnd(4, 11);
			vector<IP> q(n);
			for (bool ok = false; !ok; ) {
				for (auto& a : q) a = {rnd(-R, R), rnd(-R, R),
					slab ? rnd(0, 6) : rnd(-R, R)};
				ok = true;
				rep(i,0,n) rep(j,0,i) rep(k,0,j) rep(l,0,k)
					if (det(q[i], q[j], q[k], q[l]) == 0) ok = false;
			}
			vector<P3> p;
			for (auto& a : q) p.emplace_back(T(a[0] + off), T(a[1] - off),
				T(a[2] + off));
			vector<P3> a = p, b = p;
			auto fa = hull3d(a); check("hull3d (slow)", a, fa);
			auto fb = hull3dFast(b); check("hull3dFast", b, fb);
			cases++;
		}
	cout << "Tests passed! (" << cases << " point sets)" << endl;
}
