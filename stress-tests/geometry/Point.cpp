#include "../utilities/template.h"

#include "../../content/geometry/Point.h"

typedef __int128 L;
typedef complex<double> C;
const double PI = acos(-1.0);

int main() {
	mt19937_64 rng(7);
	auto ri = [&](ll a, ll b) { return uniform_int_distribution<ll>(a, b)(rng); };
	auto rd = [&](double a, double b) { return uniform_real_distribution<double>(a, b)(rng); };
	assert(sgn(-5) == -1 && sgn(0) == 0 && sgn(7LL) == 1 && sgn(-0.5) == -1 && sgn(0.0) == 0);
	{
		Point<ll> z;
		assert(z.x == 0 && z.y == 0);
		Point<ll> o(3);
		assert(o.x == 3 && o.y == 0);
		ostringstream os;
		os << Point<ll>(-3, 4) << Point<double>(0.5, -2);
		assert(os.str() == "(-3,4)(0.5,-2)");
	}
	// Point<ll>: exact, small values (many ties) and |x|,|y| <= 1e9
	for (ll M : {2LL, 10LL, 1000000000LL}) rep(it,0,1000000) {
		typedef Point<ll> P;
		ll ax = ri(-M, M), ay = ri(-M, M), bx = ri(-M, M), by = ri(-M, M);
		ll cx = ri(-M, M), cy = ri(-M, M), d = ri(-M, M);
		P a(ax, ay), b(bx, by), c(cx, cy);
		assert((a < b) == (make_pair(ax, ay) < make_pair(bx, by)));
		assert((a == b) == (ax == bx && ay == by));
		P s = a + b, t = a - b, m = a * d;
		assert(s.x == ax + bx && s.y == ay + by);
		assert(t.x == ax - bx && t.y == ay - by);
		assert((L)m.x == (L)ax * d && (L)m.y == (L)ay * d);
		if (d) {
			P q = a / d;
			assert(q.x == ax / d && q.y == ay / d);
		}
		assert((L)a.dot(b) == (L)ax * bx + (L)ay * by);
		assert((L)a.cross(b) == (L)ax * by - (L)ay * bx);
		assert((L)a.cross(b, c) == (L)(bx - ax) * (cy - ay) - (L)(by - ay) * (cx - ax));
		assert((L)a.dist2() == (L)ax * ax + (L)ay * ay);
		assert(abs(a.dist() - hypot((double)ax, (double)ay)) <= 1e-9 * (double)M);
		P p = a.perp();
		assert(p.x == -ay && p.y == ax);
		if (ax || ay) assert(abs(a.angle() - arg(C((double)ax, (double)ay))) < 1e-12);
	}
	// Point<double> against std::complex
	for (double M : {1.0, 1e3, 1e9}) rep(it,0,1000000) {
		typedef Point<double> P;
		P a(rd(-M, M), rd(-M, M)), b(rd(-M, M), rd(-M, M));
		if (it % 10 == 0) a = P((double)ri(-2, 2), (double)ri(-2, 2));
		C ca(a.x, a.y), cb(b.x, b.y);
		double len = abs(ca), eps = 1e-9 * max(1.0, M);
		assert(abs(a.dist() - len) <= eps);
		assert(abs(a.dist2() - norm(ca)) <= eps * max(1.0, M));
		assert(abs(a.dot(b) - (conj(ca) * cb).real()) <= eps * max(1.0, M));
		assert(abs(a.cross(b) - (conj(ca) * cb).imag()) <= eps * max(1.0, M));
		if (len == 0) continue;
		double ang = a.angle();
		assert(-PI <= ang && ang <= PI && abs(ang - arg(ca)) < 1e-12);
		P u = a.unit();
		assert(abs(u.dist() - 1) < 1e-12 && abs(u.x - a.x / len) < 1e-12 && abs(u.y - a.y / len) < 1e-12);
		P n = a.normal(), pp = a.perp();
		assert(abs(n.dist() - 1) < 1e-12 && abs(n.dot(a)) <= eps && a.cross(n) > 0);
		assert(pp.x == -a.y && pp.y == a.x);
		double th = rd(-10, 10);
		P r = a.rotate(th);
		C cr = ca * polar(1.0, th);
		assert(abs(r.x - cr.real()) <= eps && abs(r.y - cr.imag()) <= eps);
	}
	{
		typedef Point<double> P;
		P r = P(1, 0).rotate(PI / 2);
		assert(abs(r.x) < 1e-15 && abs(r.y - 1) < 1e-15); // ccw
		assert(P(-1, 0).angle() == PI && P(0, -1).angle() == -PI / 2);
	}
	cout<<"Tests passed!"<<endl;
}
