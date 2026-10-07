#include <bits/stdc++.h>

#define all(x) begin(x), end(x)
typedef long long ll;
using namespace std;

#define rep(i, a, b) for (int i = a; i < (b); ++i)
#define all(x) begin(x), end(x)
#define sz(x) (int)(x).size()
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;

#include "../../content/geometry/Point.h"
#include "../../content/geometry/sideOf.h"
#include "../../content/geometry/PolygonArea.h"
#include "../../content/geometry/PolygonUnion.h"
#include "../utilities/genPolygon.h"
#include "../utilities/random.h"

namespace blackhorse {

using db = double;
const db eps = 1e-8;

struct pt {
    db x, y;
    pt(db x = 0, db y = 0) : x(x), y(y) {}
};

inline int sgn(db x) { return (x > eps) - (x < -eps); }

pt operator-(pt p1, pt p2) { return pt(p1.x - p2.x, p1.y - p2.y); }

db vect(pt p1, pt p2) { return p1.x * p2.y - p1.y * p2.x; }

db scal(pt p1, pt p2) { return p1.x * p2.x + p1.y * p2.y; }

db polygon_union(vector<pt> poly[], int n) {
    auto ratio = [](pt A, pt B, pt O) {
        return !sgn(A.x - B.x) ? (O.y - A.y) / (B.y - A.y) : (O.x - A.x) / (B.x - A.x);
    };
    db ret = 0;
    for (int i = 0; i < n; ++i) {
        for (size_t v = 0; v < poly[i].size(); ++v) {
            pt A = poly[i][v], B = poly[i][(v + 1) % poly[i].size()];
            vector<pair<db, int>> segs;
            segs.emplace_back(0, 0), segs.emplace_back(1, 0);
            for (int j = 0; j < n; ++j)
                if (i != j) {
                    for (size_t u = 0; u < poly[j].size(); ++u) {
                        pt C = poly[j][u], D = poly[j][(u + 1) % poly[j].size()];
                        int sc = sgn(vect(B - A, C - A)), sd = sgn(vect(B - A, D - A));
                        if (!sc && !sd) {
                            if (sgn(scal(B - A, D - C)) > 0 && i > j) {
                                segs.emplace_back(ratio(A, B, C), 1), segs.emplace_back(ratio(A, B, D), -1);
                            }
                        } else {
                            db sa = vect(D - C, A - C), sb = vect(D - C, B - C);
                            if (sc >= 0 && sd < 0)
                                segs.emplace_back(sa / (sa - sb), 1);
                            else if (sc < 0 && sd >= 0)
                                segs.emplace_back(sa / (sa - sb), -1);
                        }
                    }
                }
            sort(segs.begin(), segs.end());
            db pre = min(max(segs[0].first, 0.0), 1.0), now, sum = 0;
            int cnt = segs[0].second;
            for (size_t j = 1; j < segs.size(); ++j) {
                now = min(max(segs[j].first, 0.0), 1.0);
                if (!cnt)
                    sum += now - pre;
                cnt += segs[j].second;
                pre = now;
            }
            ret += vect(A, B) * sum;
        }
    }
    return ret / 2;
}
} // namespace blackhorse

namespace approximate {
#include "../../content/geometry/InsidePolygon.h"
double polygonUnion(vector<vector<P>> &polygons, int lim) {
    int cnt = 0;
    int total = 0;
    for (double y = -lim + 1e-5; y < lim; y += lim / 500.0) {
        for (double x = -lim + 1.1e-5; x < lim; x += lim / 500.0) {
            total++;
            for (auto &i : polygons) {
                if (inPolygon(i, P(x, y))) {
                    cnt++;
                    break;
                }
            }
        }
    }
    return lim * lim * 4 * cnt / double(total);
}
} // namespace approximate

namespace lovelive {
#define re real
#define im imag
#define pb push_back
#define fir first
#define sec second
typedef double db;
const db pi = acos(db(-1));
inline int sgn(db x) { return (x > 1e-8) - (x < -1e-8); }

typedef complex<db> cpoi;
db polygon_union(vector<cpoi> py[], int n) {
    auto ratio = [](cpoi &a, cpoi &b, cpoi &c) {
        cpoi x = b - a, y = c - a;
        if (sgn(re(x)) == 0)
            return im(y) / im(x);
        return re(y) / re(x);
    };
    db ret = 0;
    for (int i = 0; i < n; ++i)
        for (size_t v = 0; v < py[i].size(); ++v) {
            cpoi a = py[i][v], b = py[i][(v + 1) % py[i].size()];
            vector<pair<db, int>> segs = {{0, 0}, {1, 0}};
            for (int j = 0; j < n; ++j)
                if (i != j)
                    for (size_t u = 0; u < py[j].size(); ++u) {
                        cpoi c = py[j][u], d = py[j][(u + 1) % py[j].size()];
                        int sc = sgn(im(conj(b - a) * (c - a)));
                        int sd = sgn(im(conj(b - a) * (d - a)));
                        if (!sc && !sd) {
                            if (sgn(re(conj(b - a) * (d - c))) > 0 && i > j) {
                                segs.pb({ratio(a, b, c), +1});
                                segs.pb({ratio(a, b, d), -1});
                            }
                        } else {
                            db sa = im(conj(d - c) * (a - c));
                            db sb = im(conj(d - c) * (b - c));
                            if (sc >= 0 && sd < 0)
                                segs.pb({sa / (sa - sb), 1});
                            else if (sc < 0 && sd >= 0)
                                segs.pb({sa / (sa - sb), -1});
                        }
                    }
            sort(segs.begin(), segs.end());
            db pre = min(max(segs[0].fir, 0.0), 1.0);
            db cur, sum = 0;
            int cnt = segs[0].sec;
            for (size_t j = 1; j < segs.size(); ++j) {
                cur = min(max(segs[j].fir, 0.0), 1.0);
                if (!cnt)
                    sum += cur - pre;
                cnt += segs[j].sec;
                pre = cur;
            }
            ret += im(conj(a) * b) * sum;
        }
    ret = abs(ret) * 0.5;
    return ret;
}
} // namespace lovelive

P randPt(int lim) { return P(randRange(-lim, lim), randRange(-lim, lim)); }

P rndUlp(int lim, long long ulps = 5) { return P(randNearIntUlps(lim, ulps), randNearIntUlps(lim, ulps)); }

P rndEps(int lim, double eps) { return P(randNearIntEps(lim, eps), randNearIntEps(lim, eps)); }

namespace slab {
typedef long double ld;
mt19937 rng(12345);
int ri(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }

// Independent oracle: vertical slab decomposition. Between two consecutive
// event x-coordinates (vertices and edge crossings) no edges cross, so the
// covered length is linear in x and equals its value at the slab midpoint.
ld slabUnion(const vector<vector<P>>& polys) {
	vector<array<ld,4>> es;
	vector<ld> xs;
	for (auto& p : polys) rep(i,0,sz(p)) {
		P a = p[i], b = p[(i+1)%sz(p)];
		es.push_back({a.x, a.y, b.x, b.y});
		xs.push_back(a.x);
	}
	rep(i,0,sz(es)) rep(j,0,i) {
		auto& e = es[i]; auto& f = es[j];
		ld dx1 = e[2]-e[0], dy1 = e[3]-e[1], dx2 = f[2]-f[0], dy2 = f[3]-f[1];
		ld den = dx1*dy2 - dy1*dx2;
		if (den == 0) continue;
		ld t = ((f[0]-e[0])*dy2 - (f[1]-e[1])*dx2) / den;
		if (t > 0 && t < 1) xs.push_back(e[0] + t*dx1);
	}
	sort(all(xs));
	ld area = 0;
	rep(k,0,sz(xs)-1) {
		ld x0 = xs[k], x1 = xs[k+1], xm = (x0+x1)/2;
		if (x1 - x0 < 1e-13L) continue;
		vector<pair<ld,ld>> iv;
		for (auto& p : polys) {
			vector<ld> ys;
			rep(i,0,sz(p)) {
				P a = p[i], b = p[(i+1)%sz(p)];
				if (a.x > b.x) swap(a, b);
				if (a.x < xm && xm < b.x)
					ys.push_back(a.y + (ld)(b.y-a.y) * (xm-a.x) / ((ld)b.x-a.x));
			}
			sort(all(ys));
			assert(sz(ys) % 2 == 0);
			for (int i = 0; i < sz(ys); i += 2) iv.emplace_back(ys[i], ys[i+1]);
		}
		sort(all(iv));
		ld len = 0, hi = -1e30L;
		for (auto [l, r] : iv) {
			if (l > hi) len += r - l, hi = r;
			else if (r > hi) len += r - hi, hi = r;
		}
		area += len * (x1 - x0);
	}
	return area;
}

typedef Point<ll> Q;
bool properOrTouch(Q a, Q b, Q c, Q d) { // closed segments intersect?
	auto o = [](Q p, Q q, Q r) { return sgn(p.cross(q, r)); };
	auto on = [](Q s, Q e, Q p) { return p.cross(s,e)==0 && (s-p).dot(e-p)<=0; };
	if (o(a,b,c)*o(a,b,d) < 0 && o(c,d,a)*o(c,d,b) < 0) return 1;
	return on(a,b,c) || on(a,b,d) || on(c,d,a) || on(c,d,b);
}
// simple polygon, positive area, distinct vertices. strict: no 180-degree angles
bool isSimple(const vector<Q>& p, bool strict) {
	int n = sz(p);
	if (n < 3) return 0;
	ll A = 0;
	rep(i,0,n) A += p[i].cross(p[(i+1)%n]);
	if (A <= 0) return 0;
	rep(i,0,n) rep(j,0,i) if (p[i] == p[j]) return 0;
	rep(i,0,n) {
		Q a = p[i], b = p[(i+1)%n], c = p[(i+2)%n];
		if (a.cross(b, c) == 0) {
			if (strict || (a-b).dot(c-b) > 0) return 0;
		}
		rep(j,0,i) {
			if ((j+1)%n == i || (i+1)%n == j) continue;
			if (properOrTouch(a, b, p[j], p[(j+1)%n])) return 0;
		}
	}
	return 1;
}
vector<Q> randSimple(int n, int lo, int hi, bool strict) {
	for (;;) {
		vector<Q> p(n);
		for (auto& q : p) q = Q(ri(lo,hi), ri(lo,hi));
		if (!isSimple(p, strict)) { reverse(all(p)); if (!isSimple(p, strict)) continue; }
		return p;
	}
}
vector<P> conv(const vector<Q>& p) { vector<P> r; for (auto q : p) r.emplace_back((double)q.x, (double)q.y); return r; }


void check(vector<vector<P>> polys) {
	ld want = slabUnion(polys);
	double got = polyUnion(polys);
	// rounding error is about eps * (largest coordinate) * (extent)
	ld scale = 1, lo = 1e30L, hi = -1e30L;
	for (auto& p : polys) for (auto q : p) {
		scale = max<ld>(scale, max(fabsl(q.x), fabsl(q.y)));
		lo = min<ld>(lo, min(q.x, q.y)), hi = max<ld>(hi, max(q.x, q.y));
	}
	if (!(fabsl(got - want) <= 1e-10L * scale * max<ld>(1, hi - lo))) {
		cout << setprecision(17) << "polyUnion " << got << ", expected " << (double)want << endl;
		for (auto& p : polys) { for (auto q : p) cout << q << ' '; cout << endl; }
		abort();
	}
}

// Small integer polygons: lots of shared/overlapping edges, touching
// vertices, identical polygons, collinear vertices and nested polygons.
void testSlab() {
	{ vector<vector<P>> none; assert(polyUnion(none) == 0); }
	check({{P(0,0), P(1,0), P(0,1)}});
	rep(it,0,150000) {
		int k = ri(1, 4);
		vector<vector<P>> ps;
		bool strict = it % 2;
		int off = it % 3 == 0 ? ri(-1000000, 1000000) : 0, mul = it % 5 == 0 ? ri(-3, 3) : 1;
		if (!mul) mul = 1000;
		rep(i,0,k) {
			int t = ri(0, 5);
			vector<Q> q;
			if (t == 0 && i) { ps.push_back(ps[ri(0,i-1)]); continue; }
			else if (t == 1) {
				int x = ri(0,3), y = ri(0,3), w = ri(1,3), h = ri(1,3);
				q = {Q(x,y), Q(x+w,y), Q(x+w,y+h), Q(x,y+h)};
			} else q = randSimple(ri(3, 6), 0, ri(2, 5), strict);
			for (auto& p : q) p = p * mul + Q(off, -off);
			rotate(q.begin(), q.begin() + ri(0, sz(q)-1), q.end());
			ps.push_back(conv(q));
		}
		check(ps);
	}
}
} // namespace slab

int numSimple = 0;
void testRandom(int n, int numPts = 10, int lim = 5, bool brute = false) {
    vector<vector<P>> polygons;
    for (int i = 0; i < n; i++) {
        vector<P> pts;
        int k = randIncl(3, numPts);
        for (int j = 0; j < k; j++) {
            pts.push_back(randPt(lim)); // rndEps(lim, 1e-10));
        }
        polygons.push_back(genPolygon(pts));
        if (polygonArea2(polygons.back()) < 0) {
            reverse(all(polygons.back()));
        }
    }
    auto val1 = polyUnion(polygons);
    vector<vector<blackhorse::pt>> polygons2;
    for (auto i : polygons) {
        vector<blackhorse::pt> t;
        for (auto j : i)
            t.push_back({j.x, j.y});
        polygons2.push_back(t);
    }
    vector<vector<lovelive::cpoi>> polygons3;
    for (auto i : polygons) {
        vector<lovelive::cpoi> t;
        for (auto j : i)
            t.push_back({j.x, j.y});
        polygons3.push_back(t);
    }
    auto val3 = blackhorse::polygon_union(polygons2.data(), sz(polygons2));
    auto val4 = lovelive::polygon_union(polygons3.data(), sz(polygons3));
    // genPolygon sometimes returns self-intersecting polygons, for which the
    // union is not well defined; only use the independent oracle on simple ones.
    auto val5 = val1;
    bool simple = true;
    for (auto &i : polygons) {
        vector<slab::Q> q;
        for (auto j : i) q.emplace_back((ll)j.x, (ll)j.y);
        simple &= slab::isSimple(q, false);
    }
    if (simple) val5 = (double)slab::slabUnion(polygons), numSimple++;
    if (abs(val1 - val3) > 1e-8 || abs(val1 - val4) > 1e-8 || abs(val1 - val5) > 1e-8 * lim * lim) {
        cout << "n=" << n << " numPts=" << numPts << " lim=" << lim << setprecision(15) << " polyUnion=" << val1 << " ref=" << val3 << " " << val4 << " slab=" << val5 << endl;
        rep(i, 0, n) {
            for (auto &x : polygons[i]) {
                cout << x << ' ';
            }
            cout << endl;
        }
        abort();
    }
}

int main() {
    // int s = (int)time(0);
    int s = 1;
    // cout << "seed " << s << endl;
    srand(s);
    for (int i = 0; i < 100; i++) {
        testRandom(2, 5, 5);
    }
    for (int i = 0; i < 100; i++) {
        testRandom(2, 10, 2);
    }
    for (int i = 0; i < 50; i++) {
        testRandom(5, 100, 5);
    }
    for (int i = 0; i < 200; i++) {
        testRandom(8, 40, 1000);
    }
    for (int i = 0; i < 2000; i++) {
        testRandom(3, 8, 4);
    }
    assert(numSimple > 1000);
    slab::testSlab();
    cout << "Tests passed!" << endl;
}
