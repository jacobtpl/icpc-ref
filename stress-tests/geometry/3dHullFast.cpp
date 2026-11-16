// Validation of 3dHullFast.h (hull3dFast): exact comparison with brute force
// and 3dHull.h in general position, exact structural checks on degenerate
// integer input (coplanar, collinear, duplicate points), non-integer doubles,
// and larger inputs.
#include "../utilities/template.h"

#include "../../content/geometry/Point3D.h"
#include "../../content/geometry/3dHull.h"
#include "../../content/geometry/3dHullFast.h"
#include "../../content/geometry/PolyhedronVolume.h"

typedef __int128 LL;
typedef array<LL,3> IP;
typedef array<int,3> Tri;

mt19937_64 rng(20261116);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }

// (b-a) x (c-a) . (d-a), on differences only
template<class V>
auto det(const V& a, const V& b, const V& c, const V& d) {
	decay_t<decltype(a[0])> u[3], v[3], w[3];
	rep(i,0,3) u[i] = b[i]-a[i], v[i] = c[i]-a[i], w[i] = d[i]-a[i];
	return (u[1]*v[2] - u[2]*v[1]) * w[0] + (u[2]*v[0] - u[0]*v[2]) * w[1]
		+ (u[0]*v[1] - u[1]*v[0]) * w[2];
}
Tri canon(Tri t) {
	while (t[0] > t[1] || t[0] > t[2]) rotate(t.begin(), t.begin()+1, t.end());
	return t;
}
template<class V> set<Tri> bruteHull(const vector<V>& p) {
	set<Tri> r; int n = sz(p);
	rep(i,0,n) rep(j,i+1,n) rep(k,j+1,n) {
		int pos = 0, neg = 0;
		rep(l,0,n) if (l != i && l != j && l != k)
			(det(p[i], p[j], p[k], p[l]) > 0 ? pos : neg)++;
		if (!pos) r.insert(canon({i, j, k}));
		else if (!neg) r.insert(canon({i, k, j}));
	}
	return r;
}
template<class V> void fail(const string& why, const vector<V>& p) {
	cout << why << "\nPoints:";
	if (sz(p) <= 30) for (auto& q : p) { cout << " (";
		rep(i,0,3) cout << (double)q[i] << (i < 2 ? "," : ")"); }
	cout << endl;
	abort();
}
template<class T> vector<Point3D<T>> conv(const vector<IP>& p) {
	vector<Point3D<T>> A;
	for (auto& q : p) A.emplace_back((T)q[0], (T)q[1], (T)q[2]);
	rep(i,0,sz(p)) if ((LL)A[i].x != p[i][0] || (LL)A[i].y != p[i][1]
		|| (LL)A[i].z != p[i][2]) fail("test bug: not representable", p);
	return A;
}

// General position: the face set is unique.
template<class T> void testGeneral(const vector<IP>& p) {
	auto A = conv<T>(p);
	set<Tri> want = bruteHull(p), got, slow;
	auto fs = hull3dFast(A);
	for (auto& f : fs) if (!got.insert(canon(f)).second)
		fail("hull3dFast: repeated face", p);
	if (got != want) fail("hull3dFast: wrong faces (" + to_string(sz(got))
		+ " vs " + to_string(sz(want)) + ")", p);
	for (auto& f : hull3d(A)) slow.insert(canon({f.a, f.b, f.c}));
	if (got != slow) fail("hull3dFast != hull3d", p);
	LL v = 0;
	for (auto& t : want) v += det(p[0], p[t[0]], p[t[1]], p[t[2]]);
	if (abs(signedPolyVolume(A, fs) * 6 - (double)v) > 1e-11 * (double)v)
		fail("hull3dFast: volume", p);
}

// Any integer input that is not coplanar: the result must be a closed,
// outward oriented triangulated surface with no point outside. Such a
// surface is the boundary of the hull. O(faces * n).
template<class T> void testExact(const vector<IP>& p) {
	auto fs = hull3dFast(conv<T>(p));
	int n = sz(p);
	map<pii, int> edges;
	set<int> verts;
	LL vol = 0;
	for (auto& f : fs) {
		rep(j,0,3) {
			if (f[j] < 0 || f[j] >= n) fail("hull3dFast: bad index", p);
			if (edges[{f[j], f[(j+1)%3]}]++) fail("hull3dFast: edge twice", p);
			verts.insert(f[j]);
		}
		LL nx = 0; // zero-area faces are not allowed
		rep(c,0,3) { IP e{}; e[c] = 1;
			IP a = p[f[0]], b = p[f[1]], d = p[f[2]];
			rep(x,0,3) b[x] -= a[x], d[x] -= a[x];
			LL cr = b[(c+1)%3]*d[(c+2)%3] - b[(c+2)%3]*d[(c+1)%3];
			nx |= cr != 0; }
		if (!nx) fail("hull3dFast: degenerate face", p);
		rep(i,0,n) if (det(p[f[0]], p[f[1]], p[f[2]], p[i]) > 0)
			fail("hull3dFast: point outside", p);
		vol += det(p[fs[0][0]], p[f[0]], p[f[1]], p[f[2]]);
	}
	for (auto& [e, c] : edges) if (!edges.count({e.second, e.first}))
		fail("hull3dFast: open edge", p);
	if (sz(verts) - sz(edges) / 2 + sz(fs) != 2 || vol <= 0)
		fail("hull3dFast: not a sphere", p);
}

bool coplanarAll(const vector<IP>& p) {
	int n = sz(p);
	rep(i,0,n) rep(j,0,i) rep(k,0,j) rep(l,0,k)
		if (det(p[i], p[j], p[k], p[l]) != 0) return false;
	return true;
}
bool generalPosition(const vector<IP>& p) {
	int n = sz(p);
	rep(i,0,n) rep(j,0,i) rep(k,0,j) rep(l,0,k)
		if (det(p[i], p[j], p[k], p[l]) == 0) return false;
	return true;
}

int main() {
	int cases = 0;
	const LL E12 = 1000000000000LL;
	// general position, all number types, far from the origin
	rep(it,0,400) for (ll R : {4LL, 50LL, 1000LL, 50000LL}) rep(slab,0,2) {
		int n = (int)rnd(4, 12);
		vector<IP> p(n);
		do for (auto& q : p) q = {rnd(-R, R), rnd(-R, R),
			slab ? rnd(0, 6) : rnd(-R, R)};
		while (!generalPosition(p));
		for (ll off : {0LL, 1000000000LL, 4000000000000000LL}) {
			auto q = p;
			for (auto& a : q) a[0] += off, a[1] -= off, a[2] += off;
			testGeneral<double>(q); testGeneral<long double>(q);
			testGeneral<ll>(q); testGeneral<LL>(q);
			testExact<ll>(q);
			cases++;
		}
		for (auto& a : p) rep(c,0,3) a[c] = a[c] * (E12 / R) + E12 * E12;
		testGeneral<LL>(p); testExact<LL>(p); // __int128: D <= 2e12
		cases++;
	}
	// degenerate integer input: tiny ranges, duplicates, lattice solids
	rep(it,0,30000) {
		int n = (int)rnd(4, 40), R = (int)rnd(1, 4);
		vector<IP> p(n);
		for (auto& q : p) q = {rnd(0, R), rnd(0, R), rnd(0, R)};
		if (it % 3 == 0) rep(i,0,n/2) p[rnd(0, n-1)] = p[rnd(0, n-1)];
		if (it % 5 == 0) for (auto& q : p) q[2] = q[2] > 0; // prism
		if (it % 7 == 0) for (auto& q : p) q[0] = q[1]; // then one apex
		if (it % 7 == 0) p[rnd(0, n-1)] = {R + 1, 0, 0};
		if (n <= 12 ? coplanarAll(p) : false) continue;
		if (n > 12) { // make sure a tetrahedron exists
			p[0] = {0,0,0}; p[1] = {R,0,0}; p[2] = {0,R,0}; p[3] = {0,0,R};
			shuffle(all(p), rng);
		}
		ll off = it % 2 ? 1000000000000LL : 0;
		for (auto& a : p) a[0] += off, a[1] -= off, a[2] += off;
		testExact<ll>(p); testExact<LL>(p); testExact<double>(p);
		cases++;
	}
	{ // full lattice cube and octahedron with interior and face points
		vector<IP> p;
		rep(x,0,7) rep(y,0,7) rep(z,0,7) p.push_back({x, y, z});
		rep(r,0,20) { shuffle(all(p), rng); testExact<ll>(p); }
		p.clear();
		rep(x,-6,7) rep(y,-6,7) rep(z,-6,7)
			if (abs(x) + abs(y) + abs(z) <= 6) p.push_back({x, y, z});
		rep(r,0,20) { shuffle(all(p), rng); testExact<ll>(p); }
		// lifted grid (Delaunay of a lattice): every facet is degenerate
		p.clear();
		rep(x,-10,11) rep(y,-10,11) p.push_back({x, y, x*x + y*y});
		rep(r,0,20) { shuffle(all(p), rng); testExact<ll>(p); }
		cases += 60;
	}
	// larger inputs: ball, sphere-like, lifted paraboloid
	rep(it,0,30) {
		int n = 3000; ll R = it % 3 == 2 ? 1000 : 100000;
		vector<IP> p(n);
		for (auto& q : p) {
			q = {rnd(-R, R), rnd(-R, R), rnd(-R, R)};
			if (it % 3 == 1) { // near a sphere: almost all on the hull
				double d = sqrt(double(q[0]*q[0] + q[1]*q[1] + q[2]*q[2])) + 1;
				rep(c,0,3) q[c] = (LL)((double)q[c] * (double)R / d);
			}
			if (it % 3 == 2) q[2] = q[0]*q[0] + q[1]*q[1];
		}
		if (it % 3 == 2) testExact<ll>(p); else testExact<LL>(p);
		cases++;
	}
	{ // non-integer doubles against brute force in long double
		typedef array<long double,3> DP;
		uniform_real_distribution<double> U(-1, 1);
		rep(it,0,40000) {
			int n = (int)rnd(4, 12), kind = it % 5;
			double off = kind == 1 ? 1e6 : 0, sc = kind == 2 ? 1e-20 :
				kind == 3 ? 1e20 : 1;
			vector<P3> A; vector<DP> p;
			rep(i,0,n) {
				double x = U(rng), y = U(rng), z = U(rng);
				if (kind == 4) { double d = sqrt(x*x + y*y + z*z);
					x /= d, y /= d, z /= d; }
				A.emplace_back(x*sc + off, y*sc - off, z*sc + off);
				p.push_back({A[i].x, A[i].y, A[i].z});
			}
			set<Tri> got;
			for (auto& f : hull3dFast(A)) got.insert(canon(f));
			if (got != bruteHull(p)) fail("hull3dFast<double>: faces", p);
			cases++;
		}
	}
	// degenerate input in doubles (lattice points scaled, rotated and
	// moved, so nothing is exact): with eps the result must still be a
	// closed surface with the exact volume and no point clearly outside
	rep(it,0,6000) {
		int n = (int)rnd(4, it % 10 ? 40 : 300), R = (int)rnd(1, 20);
		vector<IP> p(n);
		for (auto& q : p) q = {rnd(-R, R), rnd(-R, R), rnd(-R, R)};
		if (it % 3 == 0) for (auto& q : p) q[rnd(0, 2)] = it & 1 ? R : -R;
		if (it % 3 == 1) for (auto& q : p) q[2] = q[0]*q[0] + q[1]*q[1];
		p[0] = {0,0,0}; p[1] = {R,0,0}; p[2] = {0,R,0}; p[3] = {0,0,R};
		shuffle(all(p), rng);
		LL v6 = 0;
		for (auto& f : hull3dFast(conv<ll>(p)))
			v6 += det(p[0], p[f[0]], p[f[1]], p[f[2]]);
		double c = cos(0.7), s = sin(0.7), c2 = cos(1.9), s2 = sin(1.9);
		double sc = it % 2 ? 0.1 : 1, off = it % 4 < 2 ? 0 : 1e6;
		vector<P3> A;
		for (auto& q : p) {
			double x = (double)q[0], y = (double)q[1], z = (double)q[2];
			if (it % 8 < 4) { double u = c*x - s*y, w = s*x + c*y;
				x = u, y = c2*w - s2*z, z = s2*w + c2*z; }
			A.emplace_back(x*sc + off, y*sc - off, z*sc + off/3);
		}
		double D = 0;
		for (P3 a : A) D = max({D, abs(a.x - A[0].x), abs(a.y - A[0].y),
			abs(a.z - A[0].z)});
		auto fs = hull3dFast(A, 1e-9 * D*D*D);
		map<pii, int> edges;
		double vol = 0;
		for (auto& f : fs) {
			rep(j,0,3) if (edges[{f[j], f[(j+1)%3]}]++)
				fail("hull3dFast eps: edge twice", p);
			P3 a = A[f[0]], nr = (A[f[1]] - a).cross(A[f[2]] - a);
			if (!(nr.dist() > 0)) fail("hull3dFast eps: degenerate face", p);
			for (P3 q : A) if (nr.dot(q - a) > 1e-7 * D * nr.dist())
				fail("hull3dFast eps: point outside", p);
			vol += (a - A[0]).cross(A[f[1]] - A[0]).dot(A[f[2]] - A[0]);
		}
		for (auto& [e, k] : edges) if (!edges.count({e.second, e.first}))
			fail("hull3dFast eps: open edge", p);
		double want = (double)v6 * sc*sc*sc;
		if (abs(vol - want) > 1e-7 * want) fail("hull3dFast eps: volume", p);
		cases++;
	}
	{ // an int eps must not make the tests run in int
		vector<P3> A{P3(0,0,0), P3(.5,0,0), P3(0,.5,0), P3(0,0,.5)};
		assert(sz(hull3dFast(A, 0)) == 4);
	}
	{ // speed: 2e5 points on a sphere and 2e5 lifted grid points
		auto t0 = chrono::steady_clock::now();
		int n = 200000;
		vector<P3> A;
		normal_distribution<double> N;
		rep(i,0,n) { P3 q(N(rng), N(rng), N(rng)); A.push_back(q.unit()); }
		assert(sz(hull3dFast(A)) == 2 * n - 4);
		vector<Point3D<ll>> B;
		rep(i,0,n) { ll x = rnd(-300, 300), y = rnd(-300, 300);
			B.emplace_back(x, y, x*x + y*y); }
		assert(sz(hull3dFast(B)) > 1000);
		double s = chrono::duration<double>(
			chrono::steady_clock::now() - t0).count();
		cerr << "two hulls of 2e5 points: " << s << " s" << endl;
		assert(s < 120);
	}
	cout << "Tests passed! (" << cases << " point sets)" << endl;
}
