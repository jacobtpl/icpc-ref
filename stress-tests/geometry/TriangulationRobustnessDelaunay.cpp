// Exact (integer) validation of FastDelaunay.h, Delaunay.h (triHull) and
// DelaunayTriangulation.h (delaunay): collinear, duplicate, cocircular and
// large-coordinate inputs, edge-record ownership, many tiny calls.
#include "../utilities/template.h"

#include "../../content/geometry/Point.h"
#include "../../content/geometry/Point3D.h"
#include "../../content/geometry/Delaunay.h"
#include "../../content/geometry/DelaunayTriangulation.h"

typedef __int128 LL;
typedef array<P,3> Tri;

mt19937_64 rng(20261012);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }

ll orient(P a, P b, P c) { // |coords| <= 1e9: fits in long long
	return (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x); }
// > 0 iff d is strictly inside the circumcircle of ccw (a,b,c); independent
// of the header: plain cofactor expansion of the lifted 3x3 determinant.
LL inCircle(P a, P b, P c, P d) {
	LL m[3][3]; P q[3] = {a, b, c};
	rep(i,0,3) {
		m[i][0] = q[i].x - d.x; m[i][1] = q[i].y - d.y;
		m[i][2] = (LL)q[i].x*q[i].x + (LL)q[i].y*q[i].y
			- (LL)d.x*d.x - (LL)d.y*d.y;
	}
	return m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1])
		- m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0])
		+ m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);
}
vector<P> hullOf(vector<P> p) { // strict hull of distinct sorted points
	if (sz(p) < 3) return p;
	vector<P> h;
	rep(pass,0,2) {
		size_t start = h.size();
		for (P q : p) {
			while (h.size() >= start + 2 &&
				orient(h[h.size()-2], h.back(), q) <= 0) h.pop_back();
			h.push_back(q);
		}
		h.pop_back(); reverse(all(p));
	}
	return h;
}

// Empty string iff tris is a Delaunay triangulation of the point set ps.
// full = false: only local (edge-flip) legality, for large inputs.
string check(vector<P> ps, const vector<Tri>& tris, bool full) {
	sort(all(ps)); ps.erase(unique(all(ps)), ps.end());
	bool collinear = true;
	rep(i,2,sz(ps)) if (orient(ps[0], ps[1], ps[i])) collinear = false;
	if (sz(ps) < 3 || collinear)
		return tris.empty() ? "" : "triangles for a degenerate set";
	map<pair<P,P>, P> edge; // directed edge -> opposite vertex
	set<P> used;
	LL area = 0;
	for (auto& t : tris) {
		rep(i,0,3) if (!binary_search(all(ps), t[i])) return "not an input point";
		ll o = orient(t[0], t[1], t[2]);
		if (o <= 0) return "triangle not counter-clockwise";
		area += o;
		rep(i,0,3) {
			used.insert(t[i]);
			if (!edge.insert({{t[i], t[(i+1)%3]}, t[(i+2)%3]}).second)
				return "duplicate directed edge (overlapping triangles)";
		}
		if (full) for (P q : ps) if (inCircle(t[0], t[1], t[2], q) > 0)
			return "point strictly inside a circumcircle";
	}
	if (sz(used) != sz(ps)) return "unused point";
	vector<P> h = hullOf(ps);
	LL harea = 0;
	rep(i,0,sz(h)) harea += (LL)h[i].x*h[(i+1)%sz(h)].y
		- (LL)h[i].y*h[(i+1)%sz(h)].x;
	if (area != harea) return "triangles do not tile the convex hull";
	for (auto& [e, opp] : edge) {
		auto it = edge.find({e.second, e.first});
		if (it == edge.end()) { // must lie on the hull boundary
			for (P q : ps) if (orient(e.first, e.second, q) < 0)
				return "unmatched interior edge";
		} else if (inCircle(e.first, e.second, opp, it->second) > 0)
			return "illegal edge";
	}
	return "";
}

void fail(const string& what, const string& why, const vector<P>& ps) {
	cout << what << ": " << why << "\nPoints:";
	for (P p : ps) cout << " (" << p.x << "," << p.y << ")";
	cout << endl;
	abort();
}

vector<Tri> runFast(vector<P> ps) { // needs distinct points
	sort(all(ps)); ps.erase(unique(all(ps)), ps.end());
	vector<P> t = triangulate(ps);
	if (sz(t) % 3) fail("triangulate", "size not a multiple of 3", ps);
	// ownership: every edge record must be released when the call returns
	if (!quads.empty() || H) fail("triangulate", "edge records kept", ps);
	vector<Tri> r;
	for (int i = 0; i < sz(t); i += 3) r.push_back({t[i], t[i+1], t[i+2]});
	return r;
}
vector<Tri> runKactl(vector<P> ps) {
	vector<Tri> r;
	delaunay(ps, [&](int a, int b, int c) {
		r.push_back({ps[a], ps[b], ps[c]}); });
	return r;
}
bool generalPosition(const vector<P>& p) { // delaunay()'s contract
	int n = sz(p);
	rep(i,0,n) rep(j,0,i) rep(k,0,j) {
		ll o = orient(p[i], p[j], p[k]);
		if (o == 0) return false;
		rep(l,0,k) if (inCircle(p[i], o > 0 ? p[j] : p[k],
				o > 0 ? p[k] : p[j], p[l]) == 0) return false;
	}
	return true;
}

int nKactl = 0;
void testAll(const vector<P>& ps, bool full = true) {
	string e = check(ps, runFast(ps), full);
	if (!e.empty()) fail("FastDelaunay triangulate", e, ps);
	e = check(ps, triHull(ps), full); // duplicates allowed
	if (!e.empty()) fail("Delaunay triHull", e, ps);
	if (sz(ps) <= 12 && sz(ps) >= 3 && generalPosition(ps)) {
		nKactl++;
		e = check(ps, runKactl(ps), full);
		if (!e.empty()) fail("DelaunayTriangulation delaunay", e, ps);
	}
}

vector<P> gen(int mode, int n, ll C) {
	vector<P> p;
	auto cl = [&](ll v) { return max(-C, min(C, v)); };
	if (mode == 0) rep(i,0,n) p.emplace_back(rnd(-C, C), rnd(-C, C));
	else if (mode == 1) { // small lattice far away: collinear + cocircular
		ll ox = rnd(0,1) ? C - 4 : -C, oy = rnd(0,1) ? C - 4 : -C;
		rep(i,0,n) p.emplace_back(cl(ox + rnd(0,4)), cl(oy + rnd(0,4)));
	} else if (mode == 2) { // all on one line
		ll dx = rnd(-3, 3), dy = rnd(-3, 3), k = max<ll>(1, C / 40);
		P o(rnd(-C/2, C/2), rnd(-C/2, C/2));
		rep(i,0,n) { ll t = rnd(-10, 10) * rnd(1, k);
			p.emplace_back(o.x + dx*t, o.y + dy*t); }
	} else if (mode == 3) { // cocircular: lattice points of radius 65*s
		static const int v[4][2] = {{16,63},{25,60},{33,56},{39,52}};
		ll s = max<ll>(1, C / 200), ox = rnd(-C/2, C/2), oy = rnd(-C/2, C/2);
		rep(i,0,n) {
			int j = (int)rnd(0, 4); ll a = j == 4 ? 0 : v[j][0];
			ll b = j == 4 ? 65 : v[j][1];
			if (rnd(0,1)) swap(a, b);
			if (rnd(0,1)) a = -a;
			if (rnd(0,1)) b = -b;
			if (rnd(0,9) == 0) a = 0, b = 0; // sometimes the centre
			p.emplace_back(ox + a*s, oy + b*s);
		}
	} else if (mode == 4) { // extreme corners with tiny perturbations
		rep(i,0,n) p.emplace_back(
			rnd(0,3) ? (rnd(0,1) ? C - rnd(0,9) : -C + rnd(0,9)) : rnd(-C,C),
			rnd(0,3) ? (rnd(0,1) ? C - rnd(0,9) : -C + rnd(0,9)) : rnd(-C,C));
	} else { // tight cluster far from the origin, with duplicates
		ll ox = C - 20, oy = -C + rnd(0, 5);
		rep(i,0,n) p.emplace_back(ox + rnd(0, 20), oy + rnd(0, 20));
		rep(i,0,n/4) p.push_back(p[rnd(0, n-1)]);
	}
	return p;
}

int main() {
	auto pts = [](vector<ll> v) { vector<P> p;
		for (int i = 0; i < sz(v); i += 2) p.emplace_back(v[i], v[i+1]);
		return p; };
	// Reported failures of the old double-lifting implementations.
	vector<vector<P>> fixed = {
		pts({29862,-37408, -34138,-19408, 1363,-10907}), // triHull aborted
		pts({-10000000,-9999998, 9999999,-10000000, -9999998,9999999,
			10000000,-9999999}), // delaunay(): overlapping triangles
		pts({-999999998,999999996, -999999995,999999996,
			-999999993,999999995, -999999998,999999998}),
		pts({129705572,-562710026, -305294428,-267710026,
			-250294428,-402710026, 424705572,-127710026,
			-50294428,157289974}),
		pts({0,0, 1,1, 2,2, 5,5}), pts({0,0, 0,0, 4,0, 0,4, 5,5}),
		pts({0,0, 4,0, 0,0, 0,4, 5,5, 4,0}), pts({}), pts({0,0}),
		pts({0,0, 3,1}), pts({7,7, 7,7, 7,7}), pts({0,0, 4,0, 0,4, 4,4, 2,2}),
		pts({-1000000000,-1000000000, 1000000000,-1000000000,
			1000000000,1000000000, -1000000000,1000000000, 0,0}),
	};
	for (auto& f : fixed) testAll(f);
	{ // the two headline cases, checked directly against the contract
		assert(sz(triHull(fixed[0])) == 1 && sz(runKactl(fixed[0])) == 1);
		assert(check(fixed[1], runKactl(fixed[1]), true).empty());
		assert(sz(runKactl(fixed[1])) == 2 && sz(triHull(fixed[5])) == 2);
		// delaunay() keeps accepting integer-valued doubles
		vector<Point<double>> d; vector<Tri> r;
		for (P p : fixed[1]) d.emplace_back((double)p.x, (double)p.y);
		delaunay(d, [&](int a, int b, int c) {
			r.push_back({fixed[1][a], fixed[1][b], fixed[1][c]}); });
		assert(check(fixed[1], r, true).empty());
	}
	int cases = 0;
	for (ll C : {3LL, 10LL, 1000LL, 100000LL, 10000000LL, 1000000000LL})
		rep(mode,0,6) rep(it,0,1500) {
			int n = (int)rnd(1, it % 10 == 0 ? 40 : 12);
			testAll(gen(mode, n, C)); cases++;
		}
	assert(nKactl > 5000);
	// larger inputs: local legality + tiling
	for (ll C : {300LL, 1000000000LL}) rep(mode,0,2)
		testAll(gen(mode * 4, 20000, C), false);
	{ // many tiny calls must stay cheap and must not accumulate records
		auto t0 = chrono::steady_clock::now();
		size_t tot = 0;
		rep(it,0,200000) {
			vector<P> p = {P(it, 0), P(it + 4, 1), P(it + 1, 5),
				P(it + 5, 6), P(it + 2, 2)};
			tot += triangulate(p).size();
			assert(quads.empty() && !H);
		}
		assert(tot == 200000 * 12);
		double s = chrono::duration<double>(
			chrono::steady_clock::now() - t0).count();
		cerr << "200000 five-point triangulations: " << s << " s" << endl;
		assert(s < 5);
	}
	cout << "Tests passed! (" << cases << " random cases, " << nKactl
		<< " in general position for delaunay())" << endl;
}
