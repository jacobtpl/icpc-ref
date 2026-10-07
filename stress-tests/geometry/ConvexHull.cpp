#include "../utilities/template.h"

#include "../../content/geometry/ConvexHull.h"
#include "../utilities/bench.h"

namespace old {
pair<vi, vi> ulHull(const vector<P>& S) {
	vi Q(sz(S)), U, L;
	iota(all(Q), 0);
	sort(all(Q), [&S](int a, int b){ return S[a] < S[b]; });
	for(auto &it: Q) {
#define ADDP(C, cmp) while (sz(C) > 1 && S[C[sz(C)-2]].cross(\
	S[it], S[C.back()]) cmp 0) C.pop_back(); C.push_back(it);
		ADDP(U, <=); ADDP(L, >=);
	}
	return {U, L};
}

vi convexHull(const vector<P>& S) {
	vi u, l; tie(u, l) = ulHull(S);
	if (sz(S) <= 1) return u;
	if (S[u[0]] == S[u[1]]) return {0};
	l.insert(l.end(), u.rbegin()+1, u.rend()-1);
	return l;
}
}

typedef __int128 lll;
lll crossL(P a, P b, P c) {
	return (lll)(b.x - a.x) * (c.y - a.y) - (lll)(b.y - a.y) * (c.x - a.x);
}
// p is a hull vertex iff it is not in the closed hull of the other points,
// i.e. (Caratheodory) not on a segment or in a triangle spanned by them.
vector<P> bruteHull(vector<P> pts) {
	sort(all(pts));
	pts.erase(unique(all(pts)), pts.end());
	int n = sz(pts);
	vector<P> res;
	rep(i,0,n) {
		bool inside = 0;
		rep(j,0,n) rep(k,0,j) if (i != j && i != k) {
			P p = pts[i], a = pts[j], b = pts[k];
			if (crossL(a, b, p) == 0 && (lll)(a.x - p.x) * (b.x - p.x) + (lll)(a.y - p.y) * (b.y - p.y) < 0) inside = 1;
			rep(l,0,k) if (i != l) {
				P c = pts[l];
				lll s1 = crossL(a, b, p), s2 = crossL(b, c, p), s3 = crossL(c, a, p);
				if (crossL(a, b, c) != 0 && ((s1 >= 0 && s2 >= 0 && s3 >= 0) || (s1 <= 0 && s2 <= 0 && s3 <= 0))) inside = 1;
			}
		}
		if (!inside) res.push_back(pts[i]);
	}
	return res;
}
// Strictly convex, counter-clockwise, starting at the smallest point, every
// input point inside or on it; with brute: exactly the brute-force vertices.
void checkHull(const vector<P>& pts, bool brute = true) {
	vector<P> h = convexHull(pts);
	int m = sz(h);
	if (pts.empty()) { assert(h.empty()); return; }
	assert(m >= 1 && h[0] == *min_element(all(pts)));
	set<P> in(all(pts));
	for (P p : h) assert(in.count(p));
	if (m == 2) assert(!(h[0] == h[1]));
	if (m >= 3) rep(i,0,m) assert(crossL(h[i], h[(i+1)%m], h[(i+2)%m]) > 0);
	if (m >= 3) for (P p : pts) rep(i,0,m) assert(crossL(h[i], h[(i+1)%m], p) >= 0);
	if (m <= 2) for (P p : pts) assert(crossL(h[0], h[m-1], p) == 0);
	if (brute && sz(pts) <= 9) {
		vector<P> b = bruteHull(pts);
		sort(all(h));
		assert(h == b);
	}
}
ll bigRand() { return (ll)((((unsigned long long)rand() << 31) ^ rand()) % 2000000001ULL) - 1000000000LL; }

int main() {
    const int SZ = 1e2;
    rep(t,0,100000) {
        const int GRID=1e3;
        vector<P> pts(SZ);
        rep(i,0,SZ) pts[i] = P(rand()%GRID, rand()%GRID);
        auto res = convexHull(pts);
        auto res2 = old::convexHull(pts);
        assert(sz(res) == sz(res2));
        rep(i,0,sz(res2)) {
            assert(pts[res2[i]] == res[i]);
        }
    }
	// Brute force on tiny inputs (duplicates, collinear points, n = 0, 1, 2),
	// on shifted and on near-overflow coordinates.
	assert(convexHull({}).empty());
	for (ll K : {1LL, 1000000000LL / 3}) for (ll off : {0LL, -7LL}) rep(t,0,300000) {
		int n = rand() % 9, lim = rand() % 4 + 1;
		vector<P> pts(n);
		rep(i,0,n) pts[i] = P((rand() % (2*lim+1) - lim) * K + off, (rand() % (2*lim+1) - lim) * K - off);
		checkHull(pts);
	}
	// Coordinates up to 1e9 in absolute value.
	rep(t,0,20000) {
		int n = rand() % 60;
		vector<P> pts(n);
		rep(i,0,n) pts[i] = P(bigRand(), bigRand());
		if (t % 4 == 0) rep(i,0,n) pts[i] = P(rand() % 2 ? 1000000000 : -1000000000, bigRand());
		if (t % 4 == 1) rep(i,0,n) pts[i].y = pts[i].x; // one long line
		checkHull(pts);
	}
    cout<<"Tests passed!"<<endl;
}
