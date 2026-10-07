#include "../utilities/template.h"

#include "../../content/geometry/Angle.h"

typedef long double ld;
const ld PI = acosl(-1.0L);

// independent exact direction order: quadrant, then cross product
int quad(ll x, ll y) {
	if (x > 0 && y >= 0) return 0;
	if (x <= 0 && y > 0) return 1;
	if (x < 0 && y <= 0) return 2;
	return 3;
}
bool lessOracle(Angle a, Angle b) {
	if (a.t != b.t) return a.t < b.t;
	int qa = quad(a.x, a.y), qb = quad(b.x, b.y);
	if (qa != qb) return qa < qb;
	return (__int128)a.x * b.y - (__int128)a.y * b.x > 0;
}
// continuous angle; exact enough for small coordinates
ld theta(Angle a) {
	ld r = atan2l((ld)a.y, (ld)a.x);
	if (a.y < 0) r += 2 * PI;
	return r + 2 * PI * a.t;
}
bool eq(ld a, ld b) { return fabsl(a - b) < 1e-9L; }

mt19937_64 gen(4711);
int rnd(ll lo, ll hi) { return (int)((ll)(gen() % (unsigned long long)(hi - lo + 1)) + lo); }

int main() {
	const int C = 4;
	vector<Angle> small;
	rep(x,-C,C+1) rep(y,-C,C+1) if (x || y) rep(t,-1,2) small.push_back(Angle(x, y, t));

	// operator<, exhaustive on small coordinates
	for (Angle a : small) for (Angle b : small) {
		assert((a < b) == lessOracle(a, b));
		if (!eq(theta(a), theta(b))) assert((a < b) == (theta(a) < theta(b)));
		else assert(!(a < b) && !(b < a));
	}
	// operator<, random over the whole int range and near its limits
	rep(it,0,2000000) {
		int m = it % 4;
		auto co = [&]() -> int {
			if (m == 0) return rnd(-2000000000, 2000000000);
			if (m == 1) return rnd(-3, 3);
			if (m == 2) return rnd(0, 1) ? rnd(INT_MAX - 2, INT_MAX) : rnd(-INT_MAX, -INT_MAX + 2);
			return rnd(-100000, 100000);
		};
		Angle a(co(), co(), rnd(0, 1)), b(co(), co(), rnd(0, 1));
		if (it % 8 >= 4) { // (anti)parallel
			int k = rnd(-3, 3);
			if (abs((ll)a.x * k) <= INT_MAX && abs((ll)a.y * k) <= INT_MAX) b.x = a.x * k, b.y = a.y * k;
		}
		if ((!a.x && !a.y) || (!b.x && !b.y)) continue;
		assert((a < b) == lessOracle(a, b));
	}

	// rotations
	for (Angle a : small) {
		Angle r = a.t90();
		assert(r.x == -a.y && r.y == a.x && eq(theta(r), theta(a) + PI / 2));
		r = a.t180();
		assert(r.x == -a.x && r.y == -a.y && eq(theta(r), theta(a) + PI));
		r = a.t360();
		assert(r.x == a.x && r.y == a.y && eq(theta(r), theta(a) + 2 * PI));
		assert(a < a.t90() && a.t90() < a.t180() && a.t180() < a.t360());
	}

	for (Angle a : small) for (Angle b : small) {
		// operator-
		Angle d = a - b;
		assert(d.x == a.x - b.x && d.y == a.y - b.y && d.t == a.t);
		// angleDiff: exactly theta(b) - theta(a)
		Angle e = angleDiff(a, b);
		assert(e.x || e.y);
		assert(eq(theta(e), theta(b) - theta(a)));
		// point a + vector b: nearest representation of the new direction
		if (a.x + b.x || a.y + b.y) {
			Angle s = a + b;
			assert(s.x == a.x + b.x && s.y == a.y + b.y);
			assert(fabsl(theta(s) - theta(a)) < PI + 1e-9L);
		}
		// segmentAngles, for two points in the same turn
		if (a.t == b.t) {
			pair<Angle, Angle> sa = segmentAngles(a, b);
			ld lo = theta(sa.first), hi = theta(sa.second);
			ld dif = fabsl(theta(a) - theta(b));
			dif = min(dif, 2 * PI - dif);
			assert(!(sa.second < sa.first));
			assert(eq(hi - lo, dif));
			bool fa = sa.first.x == a.x && sa.first.y == a.y, fb = sa.first.x == b.x && sa.first.y == b.y;
			bool sA = sa.second.x == a.x && sa.second.y == a.y, sB = sa.second.x == b.x && sa.second.y == b.y;
			assert((fa && sB) || (fb && sA));
		}
	}

	// walking a + b + b + ... around the origin tracks the winding number
	rep(it,0,20000) {
		Angle a(rnd(-5, 5), rnd(-5, 5));
		if (!a.x && !a.y) continue;
		ld th = theta(a);
		rep(st,0,30) {
			Angle b(rnd(-3, 3), rnd(-3, 3));
			if (a.x + b.x == 0 && a.y + b.y == 0) continue;
			// skip steps through the origin (direction change is ambiguous)
			if ((ll)a.x * b.y - (ll)a.y * b.x == 0 && (ll)a.x * (a.x + b.x) + (ll)a.y * (a.y + b.y) < 0) continue;
			Angle s = a + b;
			ld dl = atan2l((ld)a.x * s.y - (ld)a.y * s.x, (ld)a.x * s.x + (ld)a.y * s.y);
			th += dl;
			assert(eq(theta(s), th));
			a = s;
		}
	}

	// the sweep from the Usage comment
	rep(it,0,20000) {
		int n = rnd(1, 8), R = rnd(1, 6);
		vector<Angle> w;
		rep(i,0,n) {
			Angle a(rnd(-R, R), rnd(-R, R));
			if (!a.x && !a.y) continue;
			bool dup = 0; // keep directions distinct
			for (Angle b : w) if (!(a < b) && !(b < a)) dup = 1;
			if (!dup) w.push_back(a);
		}
		sort(all(w));
		n = sz(w);
		vector<Angle> v = w;
		rep(i,0,n) v.push_back(w[i].t360());
		int j = 0;
		rep(i,0,n) {
			while (v[j] < v[i].t180()) ++j;
			int cnt = 0;
			rep(k,0,n) if ((ll)w[i].x * w[k].y - (ll)w[i].y * w[k].x > 0) cnt++;
			assert(j - i - 1 == cnt);
		}
	}

#ifdef ANGLE_BIG_DIFF
	// angleDiff multiplies coordinates in int: wrong (signed overflow) as
	// soon as 2*c^2 >= 2^31, i.e. for |x|,|y| > 32767
	{
		Angle a(40000, 40000), e = angleDiff(a, a);
		assert(e.y == 0 && e.x > 0); // angle 0 expected
	}
#endif
	{
		Angle a(32767, 32767), e = angleDiff(a, a), f = angleDiff(a, a.t90());
		assert(e.y == 0 && e.x > 0 && f.x == 0 && f.y > 0);
	}
	cout << "Tests passed!" << endl;
}
