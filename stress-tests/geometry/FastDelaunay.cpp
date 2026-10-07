#include "../utilities/template.h"

// #define TEST_PERF

#include "../../content/geometry/ConvexHull.h"
#include "../../content/geometry/PolygonArea.h"

#define P P2
#include "../../content/geometry/circumcircle.h"
#undef P

P2 top(P x) { return P2((double)x.x, (double)x.y); }

struct Bumpalloc {
	char buf[450 << 20];
	size_t bufp;
	void* alloc(size_t s) {
		assert(s < bufp);
		return (void*)&buf[bufp -= s];
	}
	Bumpalloc() { reset(); }

	template<class T> T* operator=(T&& x) {
		T* r = (T*)alloc(sizeof(T));
		new(r) T(move(x));
		return r;
	}
	void reset() { bufp = sizeof buf; }
} bumpalloc;

// When not testing perf, we don't want to leak memory
#ifndef TEST_PERF
#define new bumpalloc =
#endif
#include "../../content/geometry/FastDelaunay.h"
#ifndef TEST_PERF
#undef new
#endif

template<class A, class F>
void dela(A& v, F f) {
	auto ret = triangulate(v);
	assert(sz(ret) % 3 == 0);
	map<P, int> lut;
	rep(i,0,sz(v)) lut[v[i]] = i;
	for (int a = 0; a < sz(ret); a += 3) {
		f(lut[ret[a]], lut[ret[a+1]], lut[ret[a+2]]);
	}
}

// Exact verification with 128-bit arithmetic, valid for |coords| <= 5e8.
void checkExact(vector<P> ps) {
	sort(all(ps)); ps.erase(unique(all(ps)), ps.end());
	random_shuffle(all(ps));
	bumpalloc.reset(); H = 0;
	int n = sz(ps);
	vector<P> t = triangulate(ps);
	assert(sz(t) % 3 == 0);
	bool colinear = true;
	rep(i,2,n) if (ps[0].cross(ps[1], ps[i])) colinear = false;
	if (colinear) { assert(t.empty()); return; }
	set<P> used; lll sumar = 0;
	for (int i = 0; i < sz(t); i += 3) {
		P a = t[i], b = t[i+1], c = t[i+2];
		ll ar = a.cross(b, c);
		assert(ar > 0); // counter-clockwise
		sumar += ar;
		used.insert(a); used.insert(b); used.insert(c);
		for (P p : ps) { // no point strictly inside the circumcircle
			P u = a - p, v = b - p, w = c - p;
			lll det = (lll)u.dist2() * v.cross(w) + (lll)v.dist2() * w.cross(u)
				+ (lll)w.dist2() * u.cross(v);
			assert(det <= 0);
		}
	}
	for (P p : used) assert(count(all(ps), p));
	assert(sz(used) == n);
	vector<P> hull = convexHull(ps);
	assert(sumar == polygonArea2(hull));
}
void testExact() {
	const ll B = 500000000;
	checkExact({});
	checkExact({P(1,2)});
	checkExact({P(1,2), P(-3,4)});
	checkExact({P(0,0), P(1,1), P(2,2)});
	checkExact({P(0,0), P(1,0), P(0,1)});
	checkExact({P(0,0), P(0,1), P(1,0)});
	checkExact({P(-B,-B), P(B,-B), P(B,B), P(-B,B)});
	checkExact({P(-B,-B), P(B,-B), P(B,B), P(-B,B), P(0,0), P(1,0), P(B,0)});
	// 12 cocircular lattice points on a circle of radius 5e8, plus the center
	{
		vector<P> c; ll s = B / 5;
		for (ll x : {3*s, 4*s, -3*s, -4*s}) for (ll y : {3*s, 4*s, -3*s, -4*s})
			if (x*x + y*y == B*B) c.emplace_back(x, y);
		for (ll x : {B, -B}) c.emplace_back(x, 0), c.emplace_back(0, x);
		assert(sz(c) == 12);
		rep(mask,1,1<<12) if (mask % 7 == 0) {
			vector<P> sub;
			rep(i,0,12) if (mask >> i & 1) sub.push_back(c[i]);
			checkExact(sub);
			sub.emplace_back(0, 0);
			checkExact(sub);
		}
	}
	rep(it,0,60000) {
		int N = rand() % 40 + 1, type = it % 6;
		ll C = it % 3 == 0 ? 2 : it % 3 == 1 ? 1000 : B;
		vector<P> ps;
		rep(i,0,N) {
			ll x = rand() % (2*C+1) - C, y = rand() % (2*C+1) - C;
			if (type == 1) y = 3;                   // colinear
			if (type == 2 && i) y = x;              // colinear plus one point
			if (type == 3) y = y > 0 ? C : -C;      // two parallel lines
			if (type == 4 && C > 2) x = x / (C/2) * (C/2), y = y / (C/2) * (C/2);
			if (type == 5) x = i < N/2 ? -C : x;    // vertical line + cloud
			ps.emplace_back(x, y);
		}
		checkExact(ps);
	}
	// full grids: everything is cocircular
	rep(w,1,9) rep(h,1,9) {
		vector<P> ps;
		rep(x,0,w) rep(y,0,h) ps.emplace_back(x * (B/8) - B/2, y * (B/8) - B/2);
		checkExact(ps);
	}
	// larger inputs
	rep(it,0,6) {
		ll C = it % 2 ? B : 30;
		vector<P> ps;
		rep(i,0,1500) ps.emplace_back(rand() % (2*C+1) - C, rand() % (2*C+1) - C);
		checkExact(ps);
	}
}

int main1() {
	srand(2);
	testExact();
	feenableexcept(29);
	rep(it,0,3000000) {{
		bumpalloc.reset();
		// if (it % 200 == 0) cerr << endl;
		vector<P> ps;
		int N = rand() % 20 + 1;
		int xrange = rand() % 50 + 1;
		int yrange = rand() % 50 + 1;
		rep(i,0,N) {
			ps.emplace_back(rand() % (2*xrange) - xrange, rand() % (2*yrange) - yrange);
		}

		auto coc = [&](int i, int j, int k, int l) {
			double a = (ps[i] - ps[j]).dist();
			double b = (ps[j] - ps[k]).dist();
			double c = (ps[k] - ps[l]).dist();
			double d = (ps[l] - ps[i]).dist();
			double e = (ps[i] - ps[k]).dist();
			double f = (ps[j] - ps[l]).dist();
			double q = a*c + b*d - e*f;
			return abs(q) < 1e-4;
		};

		rep(i,0,N) rep(j,0,i) {
			// identical
			if (ps[i] == ps[j]) {  goto fail; }
		}
		if (false) rep(i,0,N) rep(j,0,i) rep(k,0,j) {
			// colinear
			if (ps[i].cross(ps[j], ps[k]) == 0) {  goto fail; }
		}
		if (false) rep(i,0,N) rep(j,0,i) rep(k,0,j) rep(l,0,k) {
			// concyclic
			if (coc(i,j,k,l) || coc(i,j,l,k) || coc(i,l,j,k) || coc(i,l,k,j)) {  goto fail; }
		}

		bool allColinear = true;
		if (N >= 3) {
			rep(i,2,N) if ((ps[i] - ps[0]).cross(ps[1] - ps[0])) allColinear = false;
		}

		auto fail = [&]() {
			cout << "Points:" << endl;
			for(auto &p: ps) {
				cout << p.x << ' ' << p.y << endl;
			}

			cout << "Triangles:" << endl;
			dela(ps, [&](int i, int j, int k) {
				cout << i << ' ' << j << ' ' << k << endl;
			});

			abort();
		};

		ll sumar = 0;
		vi used(N);
		bool any = false;
		dela(ps, [&](int i, int j, int k) {
			any = true;
			used[i] = used[j] = used[k] = 1;
			ll ar = ps[i].cross(ps[j], ps[k]);
			if (ar <= 0) fail();
			sumar += ar;
			P2 c = ccCenter(top(ps[i]), top(ps[j]), top(ps[k]));
			double ra = ccRadius(top(ps[i]), top(ps[j]), top(ps[k]));
			rep(l,0,N) {
				if ((top(ps[l]) - c).dist() < ra - 1e-5) fail();
			}
		});
		if (!allColinear) {
			rep(i,0,N) if (!used[i]) fail();
		} else {
			assert(!any);
		}

		vector<P> hull = convexHull(ps);
		ll ar2 = polygonArea2(hull);
		if (ar2 != sumar) fail();

		continue; }
fail:;
	}
	cout<<"Tests passed!"<<endl;
	// cerr << endl;
	return 0;
}

int main2() {
	vector<P> ps;
	int N = 100000;
	int xrange = 20000;
	int yrange = 20000;
	rep(i,0,N) {
		ps.emplace_back(rand() % (2*xrange) - xrange, rand() % (2*yrange) - yrange);
	}
	sort(all(ps));
	ps.erase(unique(all(ps)), ps.end());

	cout << sz(ps) << endl;
	triangulate(ps);
	return 0;
}

#ifdef TEST_PERF
int main() { return main2(); }
#else
int main() { return main1(); }
#endif
