// Exact validation of 3dHull.h (hull3d) and PolyhedronVolume.h on integer
// input in general position: translation far from the origin, the integer
// instantiations, and volumes of translated solids.
#include "../utilities/template.h"

#include "../../content/geometry/Point3D.h"
#include "../../content/geometry/3dHull.h"
#include "../../content/geometry/PolyhedronVolume.h"

typedef __int128 LL;
typedef array<LL,3> IP;
typedef array<int,3> Tri;

mt19937_64 rng(20261012);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }

// (b-a) x (c-a) . (d-a), on differences only
LL det(const IP& a, const IP& b, const IP& c, const IP& d) {
	LL u[3], v[3], w[3];
	rep(i,0,3) u[i] = b[i]-a[i], v[i] = c[i]-a[i], w[i] = d[i]-a[i];
	return (u[1]*v[2] - u[2]*v[1]) * w[0] + (u[2]*v[0] - u[0]*v[2]) * w[1]
		+ (u[0]*v[1] - u[1]*v[0]) * w[2];
}
bool generalPosition(const vector<IP>& p) {
	int n = sz(p);
	rep(i,0,n) rep(j,0,i) rep(k,0,j) rep(l,0,k)
		if (det(p[i], p[j], p[k], p[l]) == 0) return false;
	return true;
}
Tri canon(Tri t) {
	while (t[0] > t[1] || t[0] > t[2]) rotate(t.begin(), t.begin()+1, t.end());
	return t;
}
// Brute force: every supporting triangle, oriented outwards.
set<Tri> bruteHull(const vector<IP>& p) {
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
LL vol6(const vector<IP>& p, const set<Tri>& fs) {
	LL v = 0;
	for (auto& t : fs) v += det(p[0], p[t[0]], p[t[1]], p[t[2]]);
	return v;
}

void fail(const string& why, const vector<IP>& p) {
	cout << why << "\nPoints:";
	for (auto& q : p) { cout << " (";
		rep(i,0,3) cout << (double)q[i] << (i < 2 ? "," : ")"); }
	cout << endl;
	abort();
}

// Runs hull3d on Point3D<T> and compares faces and volume with brute force.
template<class T>
void testHull(const vector<IP>& p, bool exactVolume) {
	vector<Point3D<T>> A;
	for (auto& q : p) A.emplace_back((T)q[0], (T)q[1], (T)q[2]);
	rep(i,0,sz(p)) if ((LL)A[i].x != p[i][0] || (LL)A[i].y != p[i][1]
		|| (LL)A[i].z != p[i][2]) fail("test bug: not representable", p);
	set<Tri> want = bruteHull(p), got;
	auto fs = hull3d(A);
	for (auto& f : fs) if (!got.insert(canon({f.a, f.b, f.c})).second)
		fail("hull3d: repeated face", p);
	if (got != want) fail("hull3d: wrong faces (" + to_string(sz(got))
		+ " vs " + to_string(sz(want)) + ")", p);
	LL v = vol6(p, want);
	if (v <= 0) fail("test bug: volume", p);
	T v6 = signedPolyVolume6(A, fs);
	double vd = signedPolyVolume(A, fs);
	if (exactVolume && (LL)v6 != v) fail("signedPolyVolume6 not exact", p);
	if (abs(vd * 6 - (double)v) > 1e-11 * (double)v)
		fail("signedPolyVolume inaccurate", p);
	vector<Tri> arr(all(got)); // array<int,3> faces (other hull headers)
	if (signedPolyVolume6(A, arr) != v6) fail("array faces differ", p);
}

vector<IP> gen(int n, ll R, bool slab, LL ox, LL oy, LL oz) {
	for (;;) {
		vector<IP> p(n);
		for (auto& q : p) q = {rnd(-R, R), rnd(-R, R),
			slab ? rnd(0, 6) : rnd(-R, R)};
		if (!generalPosition(p)) continue;
		for (auto& q : p) q[0] += ox, q[1] += oy, q[2] += oz;
		return p;
	}
}

int main() {
	{ // review case: 9 points, extent ~2000, translated by 1e13
		const ll o = 10000000000000LL;
		ll c[9][3] = {{-654,844,5}, {-772,-365,2}, {47,133,5}, {195,-903,3},
			{-618,308,5}, {56,841,3}, {-393,-932,6}, {70,-165,6}, {752,-784,2}};
		vector<IP> p; vector<P3> A;
		for (auto& q : c) p.push_back({o+q[0], o+q[1], o+q[2]}),
			A.emplace_back(double(o+q[0]), double(o+q[1]), double(o+q[2]));
		auto fs = hull3d(A);
		assert(sz(fs) == 12);
		assert(signedPolyVolume6(A, fs) == 23135123);
		assert(vol6(p, bruteHull(p)) == 23135123);
		testHull<double>(p, true);
	}
	{ // unit tetrahedron at 1e9: the absolute-coordinate formula gave -1.7e8
		struct Fc { int a, b, c; };
		vector<Fc> fs = {{0,2,1}, {0,1,3}, {0,3,2}, {1,2,3}};
		vector<Tri> fa = {{0,2,1}, {0,1,3}, {0,3,2}, {1,2,3}};
		for (double o : {0.0, 1e9, -1e9, 1e15}) {
			vector<P3> A = {P3(o,o,o), P3(o+1,o,o), P3(o,o+1,o), P3(o,o,o+1)};
			assert(signedPolyVolume(A, fs) == 1.0 / 6);
			assert(signedPolyVolume(A, fa) == 1.0 / 6);
		}
		ll o = 1000000000;
		vector<Point3D<ll>> B = {Point3D<ll>(o,o,o), Point3D<ll>(o+1,o,o),
			Point3D<ll>(o,o+1,o), Point3D<ll>(o,o,o+1)};
		assert(signedPolyVolume6(B, fa) == 1);
		assert(signedPolyVolume(B, fs) == 1.0 / 6);
		assert(signedPolyVolume(B, vector<Fc>()) == 0);
	}
	int cases = 0;
	const LL E12 = 1000000000000LL;
	rep(it,0,400) for (ll R : {4LL, 50LL, 1000LL, 50000LL}) rep(slab,0,2) {
		int n = (int)rnd(4, 11);
		// double: exact while 6D^3 < 2^53 and coordinates < 2^53
		for (ll off : {0LL, 1000000LL, 1000000000LL, 1000000000000LL,
				10000000000000LL, 100000000000000LL, 4000000000000000LL}) {
			auto p = gen(n, R, slab, off, -off, off);
			testHull<double>(p, R <= 1000);
			testHull<long double>(p, R <= 1000);
			testHull<ll>(p, true);
			testHull<LL>(p, true);
			cases++;
		}
		// long long: D <= 1e6 anywhere in range (volume sum: D <= 2e5)
		for (ll R2 : {100000LL, 500000LL}) {
			auto p = gen(n, R2, slab, 4 * E12 * 1000000, -3 * E12 * 1000000, 0);
			if (R2 == 100000) testHull<ll>(p, true);
			else { vector<Point3D<ll>> A; set<Tri> got;
				for (auto& q : p) A.emplace_back((ll)q[0], (ll)q[1], (ll)q[2]);
				for (auto& f : hull3d(A)) got.insert(canon({f.a, f.b, f.c}));
				if (got != bruteHull(p)) fail("hull3d<ll>: wrong faces", p); }
			testHull<LL>(p, true);
			cases++;
		}
		// __int128: D <= 2e12, far away
		auto p = gen(n, 1000000000000LL, slab, E12 * E12 * 1000, -E12 * E12, 7);
		testHull<LL>(p, true);
		cases++;
	}
	{ // many tiny hulls: no fixed per-call cost
		auto t0 = chrono::steady_clock::now();
		size_t tot = 0;
		vector<IP> p0 = {{0,0,0}, {5,1,0}, {1,6,1}, {2,2,7}, {2,2,2},
			{-3,-4,-5}};
		assert(generalPosition(p0));
		testHull<double>(p0, true);
		size_t each = bruteHull(p0).size();
		rep(it,0,100000) {
			vector<P3> A = {P3(it,0,0), P3(it+5,1,0), P3(it+1,6,1),
				P3(it+2,2,7), P3(it+2,2,2), P3(it-3,-4,-5)};
			tot += hull3d(A).size();
		}
		assert(tot == 100000 * each);
		double s = chrono::duration<double>(
			chrono::steady_clock::now() - t0).count();
		cerr << "100000 six-point hulls: " << s << " s" << endl;
		assert(s < 5);
	}
	cout << "Tests passed! (" << cases << " random point sets)" << endl;
}
