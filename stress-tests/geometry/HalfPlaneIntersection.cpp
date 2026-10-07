#include "../utilities/template.h"
#define pb push_back // from content/contest/template.cpp

#include "../../content/geometry/HalfPlaneIntersection.h"

// Oracle: clip a huge box by each half-plane (Sutherland-Hodgman) in long double.
typedef long double ld;
typedef Point<ld> Q;
vector<Q> clip(const vector<Q>& poly, Q p, Q d) {
	vector<Q> res; int n = sz(poly);
	rep(i,0,n) {
		Q a = poly[i], b = poly[(i+1)%n];
		ld sa = d.cross(a-p), sb = d.cross(b-p);
		if (sa >= 0) res.push_back(a);
		if ((sa < 0 && sb > 0) || (sa > 0 && sb < 0))
			res.push_back(a + (b-a)*(sa/(sa-sb)));
	}
	return res;
}
template<class V> ld area(const V& v) {
	ld a = 0;
	rep(i,0,sz(v)) a += (ld)v[i].x*v[(i+1)%sz(v)].y - (ld)v[i].y*v[(i+1)%sz(v)].x;
	return a/2;
}
vector<Q> box(ld x0, ld y0, ld x1, ld y1) {
	return {Q(x0,y0),Q(x1,y0),Q(x1,y1),Q(x0,y1)};
}
vector<Q> clipAll(vector<Q> b, const vector<Ray>& rays) {
	for (auto& r : rays) b = clip(b, Q(r.p.x,r.p.y), Q(r.dp.x,r.dp.y));
	return b;
}
void fail(const vector<Ray>& rays, ld want, ld got) {
	cerr << setprecision(15) << "expected area " << (double)want
		<< " got " << (double)got << endl;
	for (auto& r : rays) cerr << " p=" << r.p << " dp=" << r.dp << endl;
	abort();
}
void check(const vector<Ray>& rays, ld want, bool bounds, ld tol, ld atol = 0) {
	vP res = halfPlaneIsect(rays, bounds);
	ld got = area(res);
	if (!(fabsl(got - want) <= tol * max((ld)1, want) + atol)) fail(rays, want, got);
	if (want > 1e-3) { // every returned vertex must satisfy every constraint
		ld sc = sqrtl(want) + 1;
		for (auto& q : res) {
			for (auto& r : rays) if (r.dp.cross(q - r.p) <
					-(1e-6 + (double)(tol * sc)) * r.dp.dist()) fail(rays, want, got);
		}
	}
}
int rnd(int C) { return rand() % (2*C+1) - C; }

int main() {
	srand(3);
	// explicit cases
	{
		vector<Ray> sq = {{P(0,0),P(1,0)},{P(1,0),P(0,1)},{P(1,1),P(-1,0)},{P(0,1),P(0,-1)}};
		check(sq, 1, false, 1e-9);
		check(sq, 1, true, 1e-9);
		// duplicated and redundant constraints
		vector<Ray> sq2 = sq; for (auto r : sq) sq2.pb(r), sq2.pb({r.p - r.dp.perp(), r.dp * 3});
		check(sq2, 1, false, 1e-9);
		// empty: two opposite parallel half-planes (y >= 1, y <= 0)
		check({{P(0,1),P(1,0)},{P(0,0),P(-1,0)}}, 0, false, 1e-9);
		check({{P(0,1),P(1,0)},{P(0,0),P(-1,0)}}, 0, true, 1e-9);
		// empty: three parallel
		check({{P(-1,-2),P(2,-3)},{P(1,3),P(2,-3)},{P(-3,-3),P(-2,3)}}, 0, false, 1e-9);
		// empty: general position
		check({{P(1,1),P(3,1)},{P(-3,-1),P(1,1)},{P(0,0),P(-2,1)},
				{P(-2,0),P(-3,-1)},{P(3,0),P(2,2)}}, 0, false, 1e-9);
		// single point
		check({{P(0,2),P(-3,3)},{P(1,-3),P(-1,-3)},{P(3,3),P(1,3)}}, 0, false, 1e-9);
		// empty, x <= 0 is collinear with the bounding box edge x >= 0
		check({{P(2,0),P(3,-5)},{P(-4,1),P(-5,-5)},{P(5,9),P(-2,4)},
				{P(2,-4),P(3,-2)},{P(0,9),P(0,5)}}, 0, true, 1e-9);
		// whole bounding box
		check({}, 1e18L, true, 1e-9);
		check({{P(5,5),P(1,1)}}, 5e17L, true, 1e-9);
	}
	// tiny integer inputs: lots of parallel/duplicate/degenerate/empty cases
	{
		int done = 0, pos = 0;
		while (done < 300000) {
			int n = rand()%8+1, C = rand()%2 ? 3 : 10;
			vector<Ray> rays;
			rep(i,0,n) {
				P d; do d = P(rnd(C), rnd(C)); while (d == P(0,0));
				rays.pb({P(rnd(C), rnd(C)), d});
			}
			const ld B = 1e6;
			auto o = clipAll(box(-B,-B,B,B), rays);
			bool unb = 0;
			for (auto q: o) if (max(fabsl(q.x),fabsl(q.y)) > B/2) unb = 1;
			if (unb) continue; // precondition: bounded
			ld ar = area(o);
			check(rays, ar, false, 1e-7);
			done++; pos += ar > 1e-6;
		}
		assert(pos > 10000);
	}
	// same with add_bounds (unbounded inputs allowed, region clipped to [0,1e9]^2)
	rep(it,0,200000) {
		int n = rand()%7, C = rand()%2 ? 5 : 1000;
		vector<Ray> rays;
		rep(i,0,n) {
			P d; do d = P(rnd(C), rnd(C)); while (d == P(0,0));
			rays.pb({P(rnd(C) + (rand()%2 ? C : 0), rnd(C) + (rand()%2 ? C : 0)), d});
		}
		auto o = clipAll(box(0,0,1e9,1e9), rays);
		ld ext = 0; // the oracle loses absolute precision on far-away slivers
		for (auto q : o) ext = max({ext, q.x, q.y});
		check(rays, area(o), true, 1e-7, 1e-9 * ext);
	}
	// integer inputs up to 1e5 (documented precision range), n up to 60
	rep(it,0,20000) {
		int n = rand()%60+3, C = it%2 ? 100000 : 1000;
		vector<Ray> rays;
		// lines through pairs of points around a circle => bounded, usually non-empty
		rep(i,0,n) {
			double a = rand()*1e-9*6.3, b = a + 0.2 + rand()%1000*1e-3;
			P p(round(C*cos(a)), round(C*sin(a))), q(round(C*cos(b)), round(C*sin(b)));
			if (p == q) continue;
			rays.pb({p, q-p});
		}
		rep(i,0,rand()%3) {
			P d; do d = P(rnd(C), rnd(C)); while (d == P(0,0));
			rays.pb({P(rnd(C), rnd(C)), d});
		}
		const ld B = 1e8;
		auto o = clipAll(box(-B,-B,B,B), rays);
		bool unb = 0;
		for (auto q: o) if (max(fabsl(q.x),fabsl(q.y)) > B/2) unb = 1;
		if (unb) continue;
		check(rays, area(o), false, 1e-7);
	}
	// many concurrent lines and collinear opposite half-planes (zero-width
	// regions) at coordinates up to 1e5
	rep(it,0,300000) {
		int K = rand()%4+3, n = rand()%10+3, C = it%3 ? 100000 : 1000;
		vector<P> pts; rep(i,0,K) pts.push_back(P(rnd(C), rnd(C)));
		vector<Ray> rays;
		rep(i,0,n) {
			int a = rand()%K, b = rand()%K; if (pts[a] == pts[b]) continue;
			rays.pb({rand()%2 ? pts[a] : pts[b], (pts[b]-pts[a]) * (rand()%2+1)});
		}
		const ld B = 1e8;
		auto o = clipAll(box(-B,-B,B,B), rays);
		bool unb = 0;
		for (auto q: o) if (max(fabsl(q.x),fabsl(q.y)) > B/2) unb = 1;
		if (unb) continue;
		check(rays, area(o), false, 1e-7, 1e-6*C);
	}
	cout<<"Tests passed!"<<endl;
}
