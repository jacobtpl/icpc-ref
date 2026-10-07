#include "../utilities/template.h"

#include "../../content/geometry/OnSegment.h"

typedef Point<ll> P;
typedef __int128 L;

// independent oracle: collinear (128 bit) and inside the bounding box
bool slow(P s, P e, P p) {
	L c = (L)(e.x - s.x) * (p.y - s.y) - (L)(e.y - s.y) * (p.x - s.x);
	return c == 0 && min(s.x, e.x) <= p.x && p.x <= max(s.x, e.x)
		&& min(s.y, e.y) <= p.y && p.y <= max(s.y, e.y);
}

int main() {
	mt19937_64 rng(5);
	auto ri = [&](ll a, ll b) { return uniform_int_distribution<ll>(a, b)(rng); };
	// exhaustive on a small grid, against enumeration of lattice points on the segment
	const int C = 3;
	ll cnt = 0;
	rep(sx,-C,C+1) rep(sy,-C,C+1) rep(ex,-C,C+1) rep(ey,-C,C+1) {
		P s(sx, sy), e(ex, ey);
		set<P> on;
		ll g = __gcd(abs(ex - sx), abs(ey - sy));
		if (g == 0) on.insert(s);
		else rep(k,0,g+1) on.insert(s + P((ex - sx) / g * k, (ey - sy) / g * k));
		rep(px,-C-1,C+2) rep(py,-C-1,C+2) {
			P p(px, py);
			assert(onSegment(s, e, p) == (bool)on.count(p));
			assert(slow(s, e, p) == (bool)on.count(p));
			cnt += on.count(p);
		}
	}
	assert(cnt > 0);
	// large coordinates (|x|,|y| <= 1e9), points on / next to / beyond the segment
	const ll M = 1000000000;
	ll hits = 0;
	rep(it,0,2000000) {
		ll g = ri(0, 3) ? ri(1, 1000) : ri(1, 1000000000);
		ll lim = M / g;
		ll dx = ri(-lim, lim), dy = ri(-lim, lim);
		P s(ri(max(-M, -M - dx * g), min(M, M - dx * g)), ri(max(-M, -M - dy * g), min(M, M - dy * g)));
		P e = s + P(dx * g, dy * g);
		if (ri(0, 20) == 0) e = s;
		ll k = ri(0, 4) ? ri(0, g) : ri(-2, g + 2);
		P p = s + P(dx * k, dy * k);
		if (ri(0, 2) == 0) p = p + P(ri(-1, 1), ri(-1, 1));
		if (ri(0, 30) == 0) p = ri(0, 1) ? s : e;
		p.x = max(-M, min(M, p.x)), p.y = max(-M, min(M, p.y));
		bool r = slow(s, e, p);
		hits += r;
		assert(onSegment(s, e, p) == r);
		assert(onSegment(e, s, p) == r);
	}
	assert(hits > 100000);
	// Point<double> with exactly representable values
	{
		typedef Point<double> D;
		assert(onSegment(D(0, 0), D(1, 1), D(0.5, 0.5)));
		assert(!onSegment(D(0, 0), D(1, 1), D(0.5, 0.25)));
		assert(!onSegment(D(0, 0), D(1, 1), D(2, 2)));
		assert(onSegment(D(-1.5, 2), D(-1.5, 2), D(-1.5, 2)));
	}
	cout<<"Tests passed!"<<endl;
}
