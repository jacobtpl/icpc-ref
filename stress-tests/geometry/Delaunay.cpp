#include "../utilities/template.h"

// Delaunay.h does not include its dependencies; these are what it needs.
// The benq hull headers additionally need pb/mp (contest template), a `pi`
// typedef and Point3D::operator!= (neither exists in the notebook).
#define pb push_back
#define mp make_pair
#include "../../content/geometry/Point.h"
#include "../../content/geometry/Point3D.h"
typedef pair<int, int> pi;
template<class T>
bool operator!=(const Point3D<T>& a, const Point3D<T>& b) { return !(a == b); }
#include "../../content/geometry/3dHull-template-benq.h"
#include "../../content/geometry/3dHull-fast-benq.h"
#include "../../content/geometry/Delaunay.h"
#include "../../content/geometry/ConvexHull.h"

typedef __int128 lll;

// > 0 iff d is strictly inside the circumcircle of ccw triangle a, b, c.
lll inCircle(P a, P b, P c, P d) {
	a = a - d, b = b - d, c = c - d;
	return (lll)a.dist2() * b.cross(c) + (lll)b.dist2() * c.cross(a) + (lll)c.dist2() * a.cross(b);
}

// Returns an empty string iff tris is a valid Delaunay triangulation of ps:
// ccw non-degenerate triangles on input points with empty circumcircles that
// tile the convex hull and use every point.
string check(const vector<P>& ps, const vector<array<P, 3>>& tris) {
	set<P> in(all(ps)), used;
	vector<P> hull = convexHull(ps);
	lll area = 0, sum = 0;
	rep(i,0,sz(hull)) area += hull[i].cross(hull[(i + 1) % sz(hull)]);
	for (auto& t : tris) {
		rep(i,0,3) {
			if (!in.count(t[i])) return "triangle vertex is not an input point";
			used.insert(t[i]);
		}
		ll ar = t[0].cross(t[1], t[2]);
		if (ar <= 0) return "triangle is degenerate or not ccw";
		sum += ar;
		for (P q : ps) if (inCircle(t[0], t[1], t[2], q) > 0)
			return "circumcircle contains an input point";
	}
	if (sum != area) return "triangles do not cover the convex hull exactly";
	if (area > 0 && used != in) return "an input point is in no triangle";
	return "";
}

vector<P> cur;
string curWhat;
void dump(const string& err) {
	cerr << curWhat << ": " << err << "; n=" << sz(cur) << " points:";
	rep(i,0,min(sz(cur), 12)) cerr << " (" << cur[i].x << "," << cur[i].y << ")";
	cerr << endl;
}
void onAbort(int s) { dump(s == SIGSEGV ? "segfault" : "aborted"); _Exit(1); }

map<string, int> failures;
void run(vector<P> ps, const string& what) {
	cur = ps, curWhat = what;
	string err = check(ps, triHull(ps));
	if (!err.empty() && failures[what + ": " + err]++ == 0) dump(err);
}

mt19937 gen(12345);
ll rnd(ll lo, ll hi) { return lo + (ll)(gen() % (unsigned long long)(hi - lo + 1)); }
P rndP(ll lim) { return P(rnd(-lim, lim), rnd(-lim, lim)); }
vector<P> distinct(vector<P> ps) {
	sort(all(ps));
	ps.erase(unique(all(ps)), ps.end());
	shuffle(all(ps), gen);
	return ps;
}
bool collinear(const vector<P>& ps) {
	for (P p : ps) for (P q : ps) if (ps[0].cross(p, q) != 0) return 0;
	return 1;
}
const vector<P> circle = {P(5,0), P(4,3), P(3,4), P(0,5), P(-3,4), P(-4,3),
	P(-5,0), P(-4,-3), P(-3,-4), P(0,-5), P(3,-4), P(4,-3)};

// Three families that stress the predicates, with |coordinates| <= 10*K.
void scaled(ll K, int iters) {
	string name = "scale " + to_string(K);
	rep(it,0,iters) { // points on one circle
		vector<P> ps;
		P off = rndP(K);
		for (P p : circle) if (gen() % 2) ps.push_back(p * K + off);
		if (sz(ps) >= 3) run(distinct(ps), name + " concyclic");
	}
	rep(it,0,iters) { // coarse grid: many concyclic and collinear subsets
		vector<P> ps;
		rep(i,0,rnd(3, 11)) ps.push_back(rndP(3) * K);
		ps = distinct(ps);
		if (sz(ps) >= 3 && !collinear(ps)) run(ps, name + " grid");
	}
	rep(it,0,iters) { // thin triangles
		vector<P> ps;
		P d(rnd(1, K), rnd(1, K));
		rep(i,0,rnd(3, 11)) ps.push_back(d * rnd(-9, 9) + rndP(1));
		ps = distinct(ps);
		if (sz(ps) >= 3 && !collinear(ps)) run(ps, name + " near-collinear");
	}
}

int main() {
	signal(SIGABRT, onAbort), signal(SIGSEGV, onAbort);
	// Fewer than three points: nothing to triangulate.
	run({}, "n=0");
	run({P(1, 2)}, "n=1");
	run({P(1, 2), P(3, 5)}, "n=2");
	// General position: no three collinear, no four concyclic.
	rep(it,0,20000) {
		int n = (int)rnd(3, 14);
		vector<P> ps;
		rep(i,0,n) ps.push_back(rndP(1000));
		bool bad = 0;
		rep(i,0,n) rep(j,0,i) rep(k,0,j) {
			if (ps[i].cross(ps[j], ps[k]) == 0) bad = 1;
			else rep(l,0,k) if (inCircle(ps[i], ps[j], ps[k], ps[l]) == 0) bad = 1;
		}
		if (!bad) run(ps, "general position");
	}
	// Distinct points on a tiny grid: lots of collinear and concyclic subsets.
	rep(it,0,100000) {
		vector<P> ps;
		ll lim = rnd(1, 4);
		rep(i,0,rnd(3, 11)) ps.push_back(rndP(lim));
		ps = distinct(ps);
		if (sz(ps) >= 3 && !collinear(ps)) run(ps, "small grid");
	}
	// Coordinates up to 3000, the range in which the double predicates are exact.
	for (ll K : {1, 10, 100, 300}) scaled(K, 5000);
	rep(it,0,300) {
		vector<P> ps;
		rep(i,0,rnd(3, 300)) ps.push_back(rndP(rnd(1, 3000)));
		ps = distinct(ps);
		if (sz(ps) >= 3 && !collinear(ps)) run(ps, "random");
	}
	// All points on one line: there are no triangles.
	rep(it,0,20000) {
		vector<P> ps;
		P a = rndP(3), d = rndP(2);
		if (d == P(0, 0)) d = P(1, 0);
		rep(i,0,rnd(3, 8)) ps.push_back(a + d * rnd(-10, 10));
		if (it % 2) ps = distinct(ps);
		if (sz(ps) >= 3) run(ps, "collinear");
	}
	// Repeated points.
	rep(it,0,100000) {
		vector<P> ps;
		ll lim = rnd(0, 3);
		rep(i,0,rnd(3, 10)) ps.push_back(rndP(lim));
		rep(i,0,rnd(1, 4)) ps.push_back(ps[gen() % sz(ps)]);
		shuffle(all(ps), gen);
		if (!collinear(ps)) run(ps, "duplicates");
	}
#ifdef DELAUNAY_BIG
	// Opt-in: hull3dFast works on Point3D<double>, so with |coordinates| above
	// roughly 3000 the lifted orientation test (~48 C^4) is no longer exact and
	// triHull aborts in prep() or returns a wrong triangulation.
	for (ll K : {1000, 10000, 100000, 1000000, 10000000, 100000000}) scaled(K, 1000);
#endif
	for (auto& f : failures) cerr << f.second << " x " << f.first << endl;
	if (!failures.empty()) return 1;
	cout<<"Tests passed!"<<endl;
}
