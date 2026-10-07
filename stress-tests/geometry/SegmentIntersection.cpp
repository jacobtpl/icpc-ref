#include "../utilities/template.h"

#include "../../content/geometry/SegmentIntersection.h"

namespace oldImpl {
template<class P>
int segmentIntersection(const P& s1, const P& e1,
        const P& s2, const P& e2, P& r1, P& r2) {
    if (e1==s1) {
        if (e2==s2) {
            if (e1==e2) { r1 = e1; return 1; } //all equal
            else return 0; //different point segments
        } else return segmentIntersection(s2,e2,s1,e1,r1,r2);//swap
    }
    //segment directions and separation
    P v1 = e1-s1, v2 = e2-s2, d = s2-s1;
    auto a = v1.cross(v2), a1 = v1.cross(d), a2 = v2.cross(d);
    if (a == 0) { //if parallel
        auto b1=s1.dot(v1), c1=e1.dot(v1),
             b2=s2.dot(v1), c2=e2.dot(v1);
        if (a1 || a2 || max(b1,min(b2,c2))>min(c1,max(b2,c2)))
            return 0;
        r1 = min(b2,c2)<b1 ? s1 : (b2<c2 ? s2 : e2);
        r2 = max(b2,c2)>c1 ? e1 : (b2>c2 ? s2 : e2);
        return 2-(r1==r2);
    }
    if (a < 0) { a = -a; a1 = -a1; a2 = -a2; }
    if (0<a1 || a<-a1 || 0<a2 || a<-a2)
        return 0;
    r1 = s1-v1*a2/a;
    return 1;
}
}
typedef Point<double> P;
bool eq(P a, P b) {
    return (a-b).dist()<1e-8;
}
// Exact oracle for integer segments, independent of segInter: the number of
// intersection points (0, 1, or 2 = infinitely many), the two endpoints of
// the common part when it is a segment, and the single point as a fraction.
typedef Point<ll> Q;
typedef __int128 lll;
struct Res { int cnt; Q lo, hi; lll nx, ny, den; };
Res exactInter(Q a, Q b, Q c, Q d) {
	Res r{0, Q(), Q(), 0, 0, 1};
	lll den = (lll)(b.x-a.x) * (d.y-c.y) - (lll)(b.y-a.y) * (d.x-c.x);
	if (den != 0) { // a + t(b-a) = c + u(d-c)
		lll tn = (lll)(c.x-a.x) * (d.y-c.y) - (lll)(c.y-a.y) * (d.x-c.x);
		lll un = (lll)(c.x-a.x) * (b.y-a.y) - (lll)(c.y-a.y) * (b.x-a.x);
		if (den < 0) den = -den, tn = -tn, un = -un;
		if (tn < 0 || tn > den || un < 0 || un > den) return r;
		r.cnt = 1;
		r.nx = (lll)a.x * den + tn * (b.x-a.x);
		r.ny = (lll)a.y * den + tn * (b.y-a.y);
		r.den = den;
		return r;
	}
	// parallel or degenerate: all four points must be on one line
	Q p = a, q = a == b ? (c == d ? a : c) : b;
	if (a == b && !(c == d)) p = d;
	for (Q x : {a, b, c, d})
		if ((lll)(q.x-p.x) * (x.y-p.y) - (lll)(q.y-p.y) * (x.x-p.x) != 0) return r;
	Q lo = max(min(a, b), min(c, d)), hi = min(max(a, b), max(c, d));
	if (hi < lo) return r;
	r.cnt = lo == hi ? 1 : 2;
	r.lo = lo, r.hi = hi;
	r.nx = lo.x, r.ny = lo.y;
	return r;
}

mt19937_64 rng(17);
ll ri(ll a, ll b) { return uniform_int_distribution<ll>(a, b)(rng); }

void checkExact(Q a, Q b, Q c, Q d) {
	Res e = exactInter(a, b, c, d);
	typedef Point<double> PD;
	auto dbl = [](Q p) { return PD((double)p.x, (double)p.y); };
	auto ri_ = segInter(a, b, c, d);
	auto rd = segInter(dbl(a), dbl(b), dbl(c), dbl(d));
	if (sz(ri_) != e.cnt || sz(rd) != e.cnt) {
		cout << a << ' ' << b << ' ' << c << ' ' << d << " expected " << e.cnt
			<< " points, got " << sz(ri_) << " (ll) " << sz(rd) << " (double)" << endl;
		abort();
	}
	if (e.cnt == 2) {
		assert(ri_[0] == e.lo && ri_[1] == e.hi);
		assert(rd[0] == dbl(e.lo) && rd[1] == dbl(e.hi));
	} else if (e.cnt == 1) {
		long double x = (long double)e.nx / (long double)e.den, y = (long double)e.ny / (long double)e.den;
		long double sc = 1;
		for (Q p : {a, b, c, d}) sc = max<long double>(sc, (long double)max(abs(p.x), abs(p.y)));
		assert(fabsl(rd[0].x - x) <= 1e-9L * sc && fabsl(rd[0].y - y) <= 1e-9L * sc);
		// Point<ll> is only exact when the intersection is a lattice point
		if (e.nx % e.den == 0 && e.ny % e.den == 0)
			assert(ri_[0] == Q((ll)(e.nx / e.den), (ll)(e.ny / e.den)));
	}
}

void testExact() {
	// every pair of segments on a 4x4 grid
	const int G = 4;
	rep(m,0,1<<16) {
		int v[8];
		rep(i,0,8) v[i] = (m >> (2*i)) & 3;
		checkExact(Q(v[0],v[1]), Q(v[2],v[3]), Q(v[4],v[5]), Q(v[6],v[7]));
	}
	(void)G;
	// negative and large coordinates (products of three coordinates must fit)
	for (ll lim : {3LL, 50LL, 500000LL}) rep(it,0,700000) {
		Q a(ri(-lim,lim), ri(-lim,lim)), b(ri(-lim,lim), ri(-lim,lim));
		Q c(ri(-lim,lim), ri(-lim,lim)), d(ri(-lim,lim), ri(-lim,lim));
		int t = (int)ri(0, 9);
		if (t == 0) { // collinear
			Q dir = b - a;
			ll g = max<ll>(1, __gcd(abs(dir.x), abs(dir.y)));
			dir = dir / g;
			c = a + dir * ri(-g - 2, 2 * g + 2), d = a + dir * ri(-g - 2, 2 * g + 2);
		} else if (t == 1) c = ri(0, 1) ? a : b; // shared endpoint
		else if (t == 2) { // segments crossing in a lattice point
			Q x(ri(-lim/4, lim/4), ri(-lim/4, lim/4));
			Q u(ri(-lim/8, lim/8), ri(-lim/8, lim/8)), w(ri(-lim/8, lim/8), ri(-lim/8, lim/8));
			a = x - u * ri(0, 3), b = x + u * ri(0, 3);
			c = x - w * ri(0, 3), d = x + w * ri(0, 3);
		} else if (t == 3) b = a; // point segment
		else if (t == 4) d = c + (b - a); // parallel
		checkExact(a, b, c, d);
	}
}

int main() {
	testExact();
    rep(t,0,1000000) {
        const int GRID=6;
        P a(rand()%GRID, rand()%GRID), b(rand()%GRID, rand()%GRID), c(rand()%GRID, rand()%GRID), d(rand()%GRID, rand()%GRID);
        P tmp1, tmp2;
        auto res = oldImpl::segmentIntersection(a,b,c,d, tmp1, tmp2);
        auto res2 = segInter(a,b,c,d);
        if (res != sz(res2)) {
            cout<<a<<' '<<b<<' '<<c<<' '<<d<<endl;
            cout<<"old: "<<res<<" new: "<<sz(res2)<<endl;
        }
        assert(res==sz(res2));
        if (res==1) {
            assert(eq(*res2.begin(), tmp1));
        } else if (res==2) {
            vector<P> a(res2.begin(), res2.end());
            vector<P> b({tmp1, tmp2});
            sort(all(b));
            assert(eq(a[0], b[0]) && eq(a[1],b[1]));
        }
    }
    cout<<"Tests passed!"<<endl;
}
