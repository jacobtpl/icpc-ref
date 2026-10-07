#include "../utilities/template.h"

#include "../../content/geometry/PolygonArea.h"
#include "../../content/geometry/PolygonCut.h"
#include "../../content/geometry/sideOf.h"
#include "../../content/geometry/InsidePolygon.h"
#include "../../content/geometry/SegmentIntersection.h"

typedef Point<double> P;
typedef long double ld;
mt19937 rng(7);
int ri(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }
double area(const vector<P>& p) {
	double a = 0;
	rep(i,0,sz(p)) a += p[i].cross(p[(i+1)%sz(p)]);
	return a / 2;
}

// Oracle: signed area of (triangle a,b,c) intersected with the kept half-plane
// (strictly right of s->e), by case analysis on the number of kept vertices.
ld triCut(P a, P b, P c, P s, P e) {
	ld A = (ld)a.cross(b, c) / 2;
	ld d[3] = {(ld)s.cross(e, a), (ld)s.cross(e, b), (ld)s.cross(e, c)};
	int in = (d[0] < 0) + (d[1] < 0) + (d[2] < 0);
	if (in == 0) return 0;
	if (in == 3) return A;
	rep(i,0,3) { // corner i is the odd one out
		ld x = d[i], y = d[(i+1)%3], z = d[(i+2)%3];
		if ((x < 0) == (y < 0) || (x < 0) == (z < 0)) continue;
		ld corner = A * (x / (x - y)) * (x / (x - z));
		return in == 1 ? corner : A - corner;
	}
	assert(0); return 0;
}
// Works for any (even self-intersecting) polygon since both sides are linear
// in the winding number.
ld cutOracle(const vector<P>& p, P s, P e) {
	ld r = 0;
	rep(i,1,sz(p)-1) r += triCut(p[0], p[i], p[i+1], s, e);
	return r;
}

void testOracle() {
	rep(it,0,300000) {
		int n = ri(3, 8), lim = it % 3 == 0 ? 1000000 : ri(1, 6);
		vector<P> ps;
		rep(i,0,n) ps.emplace_back(ri(-lim, lim), ri(-lim, lim));
		if (it % 5 == 0) ps[ri(0, n-1)] = ps[ri(0, n-1)]; // duplicate vertex
		P s(ri(-lim, lim), ri(-lim, lim)), e(ri(-lim, lim), ri(-lim, lim));
		if (it % 7 == 0) s = ps[0], e = ps[1]; // cut along an edge
		if (it % 11 == 0) s = ps[ri(0, n-1)]; // cut through a vertex
		auto res = polygonCut(ps, s, e);
		ld want = cutOracle(ps, s, e);
		assert(fabsl(area(res) - want) <= 1e-9L * lim * lim);
		// every vertex of the result is on the kept side or on the line
		for (P q : res) assert(s.cross(e, q) <= 1e-6 * lim * lim);
		if (s == e) assert(res.empty());
	}
	// degenerate sizes
	assert(polygonCut({}, P(0,0), P(1,0)).empty());
	assert(sz(polygonCut({P(0,-1)}, P(0,0), P(1,0))) == 1);
	assert(polygonCut({P(0,1)}, P(0,0), P(1,0)).empty());
	assert(sz(polygonCut({P(0,-1), P(1,-1)}, P(0,0), P(1,0))) == 2);
	// usage example from the header: keeps y < 0
	vector<P> sq = {P(-1,-1), P(1,-1), P(1,1), P(-1,1)};
	sq = polygonCut(sq, P(0,0), P(1,0));
	assert(sz(sq) == 4 && abs(area(sq) - 2) < 1e-12);
}

// Cutting with a line that (up to rounding) contains an edge of the polygon,
// e.g. the same half-plane twice, must not change the polygon.
void testRepeatedCut() {
	{
		vector<P> p = {P(-1,-1), P(1,-1), P(1,1), P(-1,1)};
		P s(-0.1,-0.2), e(-0.7,-0.4);
		p = polygonCut(p, s, e);
		double a = area(p);
		assert(abs(a - 7 / 3.0) < 1e-12);
		p = polygonCut(p, s, e);
		assert(abs(area(p) - a) < 1e-9);
	}
	rep(it,0,500000) {
		int D = ri(0, 1) ? 1 : 10;
		auto rp = [&]() { return P(ri(-10,10) / (double)D, ri(-10,10) / (double)D); };
		vector<P> p = {P(-1,-1), P(1,-1), P(1,1), P(-1,1)};
		if (ri(0, 1)) p = {rp(), rp(), rp()};
		if (area(p) <= 0) continue;
		rep(c,0,3) {
			P s = rp(), e = rp();
			if (s == e) continue;
			ld want = cutOracle(p, s, e);
			p = polygonCut(p, s, e);
			double a0 = area(p);
			assert(fabsl(a0 - want) < 1e-9);
			P s2 = s, e2 = e; // the same line, possibly written differently
			int t = ri(0, 2);
			if (t == 1) e2 = s + (e - s) * 2;
			if (t == 2) s2 = e + (e - s), e2 = s2 + (e - s);
			p = polygonCut(p, s2, e2);
			if (abs(area(p) - a0) > 1e-9) {
				cout << setprecision(17) << "repeated cut " << s2 << ' ' << e2
					<< " changed area " << a0 << " -> " << area(p) << endl;
				abort();
			}
			for (P q : p) assert(abs(q.x) < 10.001 && abs(q.y) < 10.001);
		}
	}
}

int main() {
	testOracle();
	testRepeatedCut();
	rep(it,0,500) {
		int N = rand() % 10 + 3;
		vector<P> ps;
		rep(i,0,N) ps.emplace_back(rand() % 10 - 5, rand() % 10 - 5);
		P p(rand() % 10 - 5), q(rand() % 10 - 5);
		rep(i,0,N) rep(j,i+1,N) {
			P a = ps[i], b = ps[(i+1)%N];
			P c = ps[j], d = ps[(j+1)%N];
			P r1, r2;
			auto r = segInter(a, b, c, d);
			if (sz(r) == 2) goto fail;
			if (sz(r) == 1) {
				if (i+1 == j || (j+1) % N == i) ;
				else goto fail;
			}
		}
		if (p == q) { fail: continue; }

		int count = 0;
		const int ITS = 400000;
		rep(it,0,ITS) {
			double x = rand() / (RAND_MAX + 1.0) * 10 - 5;
			double y = rand() / (RAND_MAX + 1.0) * 10 - 5;
			if (!inPolygon(ps, P{x,y}, true)) continue;
			if (sideOf(p, q, P{x,y}) > 0) continue;
			count++;
		}
		double approxArea = (double)count / ITS * 100;

		ps = polygonCut(ps, p, q);
		double realArea = ps.empty() ? 0.0 : abs(polygonArea2(ps) / 2.0);

		// cout << setprecision(2) << fixed;
		assert(realArea - approxArea < 2e-1);
		// cout << N << ' ' << realArea << '\t' << approxArea << '\t' << realArea - approxArea << endl;

		// cerr << N << endl;
		// for(auto &x: ps) {
			// cout << x.x << ' ' << x.y << endl;
		// }
	}
	cout<<"Tests passed!"<<endl;
}
