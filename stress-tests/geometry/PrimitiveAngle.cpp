#include "../utilities/template.h"

#include "../../content/geometry/Angle.h"

typedef __int128 LL;

// independent order: full turns, then half plane, then cross
bool lessO(Angle a, Angle b) {
	int ha = a.y < 0 || (a.y == 0 && a.x < 0);
	int hb = b.y < 0 || (b.y == 0 && b.x < 0);
	if (a.t != b.t) return a.t < b.t;
	if (ha != hb) return ha < hb;
	return (LL)a.x * b.y - (LL)a.y * b.x > 0;
}

void checkDiff(Angle a, Angle b) {
	Angle d = angleDiff(a, b);
	LL ex = (LL)a.x * b.x + (LL)a.y * b.y;
	LL ey = (LL)a.x * b.y - (LL)a.y * b.x;
	assert(d.x == ex && d.y == ey);
	// position of b's direction relative to a's within one turn
	Angle a0(a.x, a.y), b0(b.x, b.y);
	assert(d.t == b.t - a.t - lessO(b0, a0));
}

int main() {
	{ // reviewed reproducer: int products overflowed
		Angle a(50000, 1), b(1, 50000);
		Angle d = angleDiff(a, b);
		assert(d.x == 100000 && d.y == 2499999999LL && d.t == 0);
		assert(d < Angle(-1, 0)); // less than 180 degrees
	}
	{ // top of the documented range, 2C^2 = 8e18
		ll C = 2000000000;
		for (ll x : {-C, C}) for (ll y : {-C, C})
			for (ll z : {-C, C}) for (ll w : {-C, C})
				checkDiff(Angle(x, y), Angle(z, w, 1));
		Angle d = angleDiff(Angle(C, C), Angle(-C, C));
		assert(d.x == 0 && d.y == 2*C*C && d.t == 0);
		Angle e = angleDiff(Angle(C, C-1), Angle(C-1, C));
		assert(e.x == 2*C*(C-1) && e.y == 2*C-1);
		// differences of such magnitude still compare exactly
		assert(e < d && !(d < e) && e < d.t180());
		Angle f(e.x - 1, e.y);
		assert(e < f && !(f < e)); // cross = e.y > 0
	}
	// exhaustive small grid against atan2 and the exact oracle
	const double PI = acos(-1);
	vector<Angle> g;
	rep(x,-4,5) rep(y,-4,5) rep(t,-1,2) if (x || y)
		g.push_back(Angle(x, y, t));
	for (Angle a : g) for (Angle b : g) {
		assert((a < b) == lessO(a, b));
		checkDiff(a, b);
		Angle d = angleDiff(a, b);
		auto full = [&](Angle v) {
			double r = atan2((double)v.y, (double)v.x);
			if (r < -1e-12) r += 2*PI;
			return r + 2*PI*v.t;
		};
		if (d.x || d.y)
			assert(abs(full(d) - (full(b) - full(a))) < 1e-9);
	}
	// random large coordinates
	mt19937_64 rng(2024);
	for (ll C : {46341LL, 1000000LL, 1000000000LL, 2000000000LL})
		rep(it,0,300000) {
			auto r = [&]() {
				return uniform_int_distribution<ll>(-C, C)(rng); };
			Angle a(r(), r(), int(rng() % 3)), b(r(), r(), int(rng() % 3));
			if (it % 4 == 0) b.x = a.x + (ll)(rng() % 3) - 1;
			if (it % 8 == 0) b.y = a.y + (ll)(rng() % 3) - 1;
			if (it % 16 == 0) b.x = -a.x, b.y = -a.y;
			if (!(a.x || a.y) || !(b.x || b.y)) continue;
			assert((a < b) == lessO(a, b));
			checkDiff(a, b);
			Angle d = angleDiff(a, b), e = angleDiff(b, a);
			if (d.x || d.y) // differences are comparable too
				assert((d < e) == lessO(d, e));
			assert(a.t180().x == -a.x && a.t90().x == -a.y);
			if (a.x + b.x == 0 && a.y + b.y == 0) continue;
			Angle s = a + Angle(b.x, b.y);
			assert(s.x == a.x + b.x && abs(s.t - a.t) <= 1);
		}
	// operator< beyond int: all long long coordinates
	ll H = 4000000000000000000LL;
	assert(Angle(H, H - 1) < Angle(H - 1, H));
	assert(!(Angle(H - 1, H) < Angle(H, H - 1)));
	assert(Angle(H, -H) .t180() < Angle(H, -H + 1).t180());
	// many tiny sorts must stay cheap
	ll sum = 0;
	rep(it,0,100000) {
		vector<Angle> v;
		rep(i,0,4) v.push_back(Angle((ll)(rng() % 7) - 3, (ll)(rng() % 7) + 1));
		sort(all(v));
		sum += v[0].x;
	}
	assert(sum != LLONG_MAX);
	cout << "Tests passed!" << endl;
}
