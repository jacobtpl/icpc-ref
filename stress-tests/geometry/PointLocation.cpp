#include "../utilities/template.h"

#include "../../content/geometry/PointLocation.h"
#include "../../content/geometry/OnSegment.h"

typedef __int128 LL;
typedef pair<P, P> Seg;

// Can s be added to a set of segments satisfying the preconditions?
bool compatible(Seg s, Seg t) {
	auto [a, b] = s; auto [c, d] = t;
	if (a == b) return false;
	if ((a == c && b == d) || (a == d && b == c)) return true;
	ll c1 = a.cross(b, c), c2 = a.cross(b, d);
	ll c3 = c.cross(d, a), c4 = c.cross(d, b);
	if (c1 == 0 && c2 == 0) { // collinear: overlap at most a point
		if (b < a) swap(a, b);
		if (d < c) swap(c, d);
		return b < c || d < a || b == c || d == a;
	}
	if (sgn(c1) * sgn(c2) < 0 && sgn(c3) * sgn(c4) < 0) return false;
	return true;
}

// Brute force straight from the definition: among non-vertical segments
// with a <= q <= b that are below q at x = q.x, take the highest one
// (ties broken by slope, since the ray leans towards +x).
// Returns all acceptable indices.
vi brute(const vector<Seg>& segs, P q) {
	vi on;
	rep(i,0,sz(segs)) if (onSegment(segs[i].first, segs[i].second, q))
		on.push_back(i);
	if (!on.empty()) return on;
	vi best;
	LL bn = 0, bd = 1; ll bsn = 0, bsd = 1; // y = bn/bd, slope = bsn/bsd
	rep(i,0,sz(segs)) {
		auto [a, b] = segs[i];
		if (b < a) swap(a, b);
		if (a.x == b.x || q < a || b < q) continue;
		ll dx = b.x - a.x, dy = b.y - a.y;
		LL num = (LL)a.y * dx + (LL)(q.x - a.x) * dy, den = dx;
		if (num >= (LL)q.y * den) continue; // not below q
		if (best.empty()) { best = {i}; bn = num, bd = den; bsn = dy, bsd = dx; continue; }
		LL diff = num * bd - bn * den;
		if (diff == 0) {
			ll sdiff = dy * bsd - bsn * dx;
			if (sdiff == 0) { best.push_back(i); continue; }
			diff = sdiff;
		}
		if (diff > 0) { best = {i}; bn = num, bd = den; bsn = dy, bsd = dx; }
	}
	if (best.empty()) best.push_back(-1);
	return best;
}

void check(vector<Seg> segs, const vector<P>& qs) {
	vi res = pointLocation(segs, qs);
	assert(sz(res) == sz(qs));
	rep(j,0,sz(qs)) {
		vi exp = brute(segs, qs[j]);
		if (!count(all(exp), res[j])) {
			cout << "segments:" << endl;
			for (auto [a, b] : segs) cout << a << " " << b << endl;
			cout << "query " << qs[j] << " got " << res[j] << " expected one of";
			for (int x : exp) cout << " " << x;
			cout << endl;
			assert(false);
		}
	}
}

P randP(int R) { return P(rand() % (2*R+1) - R, rand() % (2*R+1) - R); }

vector<Seg> genSegs(int R, int n, int mode) {
	vector<P> pool;
	rep(i,0,R+2) pool.push_back(randP(R));
	vector<Seg> segs;
	rep(it,0,n*3) {
		if (sz(segs) >= n) break;
		P a, b;
		if (mode == 1) a = pool[rand() % sz(pool)], b = pool[rand() % sz(pool)];
		else a = randP(R), b = randP(R);
		if (mode == 2) b.x = a.x; // many verticals
		if (mode == 3) b.y = a.y; // many horizontals
		if (mode == 4 && rand() % 2) b = pool[rand() % sz(pool)];
		if (a == b) continue;
		Seg s(a, b);
		bool ok = true;
		for (Seg& t : segs) ok &= compatible(s, t);
		if (ok) segs.push_back(s);
	}
	if (!segs.empty() && rand() % 3 == 0) { // duplicates
		Seg s = segs[rand() % sz(segs)];
		if (rand() % 2) swap(s.first, s.second);
		segs.push_back(s);
	}
	return segs;
}

void testSmall() {
	rep(it,0,30000) {
		int R = rand() % 4 + 1, n = rand() % 12, mode = rand() % 5;
		vector<Seg> segs = genSegs(R, n, mode);
		vector<P> qs;
		rep(x,-R-1,R+2) rep(y,-R-1,R+2) qs.emplace_back(x, y);
		rep(i,0,5) qs.push_back(randP(R));
		check(segs, qs);

		// same configuration scaled/translated to fill +-1e9
		ll K = (ll)1e9 / (R + 1), ox = rand() % 3 - 1, oy = rand() % 3 - 1;
		auto f = [&](P p) { return P(p.x * K + ox, p.y * K + oy); };
		vector<Seg> big;
		for (auto [a, b] : segs) big.emplace_back(f(a), f(b));
		vector<P> bq;
		for (P q : qs) {
			bq.push_back(f(q));
			rep(d,0,2) bq.push_back(f(q) + P(rand() % 5 - 2, rand() % 5 - 2));
		}
		for (P& q : bq) q.x = clamp(q.x, (ll)-1e9, (ll)1e9), q.y = clamp(q.y, (ll)-1e9, (ll)1e9);
		check(big, bq);
	}
}

void testMedium() {
	rep(it,0,30) {
		int R = 30 + rand() % 20, n = 80 + rand() % 40, mode = rand() % 5;
		vector<Seg> segs = genSegs(R, n, mode);
		vector<P> qs;
		rep(i,0,2000) {
			P q = randP(R + 1);
			if (!segs.empty() && rand() % 2) {
				Seg s = segs[rand() % sz(segs)];
				q = (rand() % 2 ? s.first : s.second) + P(rand() % 3 - 1, rand() % 3 - 1);
			}
			qs.push_back(q);
		}
		check(segs, qs);
	}
}

void testEdge() {
	const ll M = 1e9;
	check({}, {});
	check({}, {P(0, 0), P(M, -M)});
	check({{P(0, 0), P(5, 0)}}, {});
	rep(i,0,2) {
		vector<Seg> segs = {{P(0, 0), P(5, 0)}};
		if (i) swap(segs[0].first, segs[0].second);
		vector<P> qs;
		rep(x,-1,7) rep(y,-1,2) qs.emplace_back(x, y);
		check(segs, qs);
		segs = {{P(0, 0), P(0, 5)}};
		if (i) swap(segs[0].first, segs[0].second);
		qs.clear();
		rep(x,-1,2) rep(y,-1,7) qs.emplace_back(x, y);
		check(segs, qs);
	}
	// bounding box and diagonals at exactly +-1e9
	vector<P> c = {P(-M, -M), P(M, -M), P(M, M), P(-M, M)};
	vector<Seg> box;
	rep(i,0,4) box.emplace_back(c[i], c[(i + 1) % 4]);
	rep(i,0,4) box.emplace_back(P(0, 0), c[i]); // half-diagonals
	vector<P> qs = {P(0, 0), P(0, 1), P(0, -1), P(1, 0), P(-1, 0)};
	for (P p : c) rep(dx,-1,2) rep(dy,-1,2) {
		P q = p + P(dx, dy);
		if (abs(q.x) <= M && abs(q.y) <= M) qs.push_back(q);
	}
	rep(i,0,200) qs.push_back(P(rand() % 3 - 1, rand() % 3 - 1) * (M / 2) + P(rand() % 5 - 2, rand() % 5 - 2));
	check(box, qs);
}

int main() {
	srand(2);
	testEdge();
	testSmall();
	testMedium();
	cout << "Tests passed!" << endl;
}
