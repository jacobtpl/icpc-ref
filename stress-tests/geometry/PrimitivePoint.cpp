#include "../utilities/template.h"

#include "../../content/geometry/Point.h"
#include "../../content/geometry/Point3D.h"

typedef long double ld;
typedef Point<ll> PL;
typedef Point<double> PD;
typedef Point<ld> PE;
typedef Point3D<ll> QL;
typedef Point3D<double> QD;
typedef Point3D<ld> QE;

template<class A, class B> void same() {
	static_assert(is_same<A, B>::value, "unexpected type"); }

void types() {
	PL a(3, 4); PD b(3, 4); PE c(3, 4); Point<int> d(3, 4);
	same<decltype(a.dist()), double>();
	same<decltype(b.dist()), double>();
	same<decltype(c.dist()), ld>();
	same<decltype(d.dist()), double>();
	same<decltype(c.angle()), ld>();
	same<decltype(a.unit()), PD>();
	same<decltype(a.normal()), PD>();
	same<decltype(a.rotate(1)), PD>();
	same<decltype(b.unit()), PD>();
	same<decltype(b.rotate(1)), PD>();
	same<decltype(c.unit()), PE>();
	same<decltype(c.rotate(1)), PE>();
	same<decltype(a.dist2()), ll>();
	same<decltype(a.cross(a)), ll>();
	same<decltype(a.perp()), PL>();
	QL e(1, 2, 3); QD f(1, 2, 3); QE g(1, 2, 3);
	same<decltype(e.dist()), double>();
	same<decltype(g.dist()), ld>();
	same<decltype(g.phi()), ld>();
	same<decltype(g.theta()), ld>();
	same<decltype(e.unit()), QD>();
	same<decltype(e.normal(e)), QD>();
	same<decltype(e.rotate(1, e)), QD>();
	same<decltype(f.unit()), QD>();
	same<decltype(f.rotate(1, f)), QD>();
	same<decltype(g.rotate(1, g)), QE>();
	same<decltype(e.cross(e)), QL>();
}

// exact user types must still instantiate (no sqrt needed)
struct Frac {
	ll n, d;
	Frac(ll n = 0, ll d = 1) : n(n), d(d) {}
	Frac operator+(Frac o) const { return {n*o.d + o.n*d, d*o.d}; }
	Frac operator-(Frac o) const { return {n*o.d - o.n*d, d*o.d}; }
	Frac operator*(Frac o) const { return {n*o.n, d*o.d}; }
	bool operator==(Frac o) const { return n*o.d == o.n*d; }
	bool operator<(Frac o) const { return n*o.d < o.n*d; }
};

void operators() {
	Point<Frac> fa(Frac(1, 2), Frac(1, 3)), fb(Frac(1, 2), Frac(2, 3));
	assert((fa + fb).x == Frac(1) && fa.cross(fb) == Frac(1, 6));
	Point<__int128> wa(1, 2), wb(3, 5);
	assert(wa.cross(wb) == -1 && (wa += wb) == Point<__int128>(4, 7));
	Point3D<__int128> wc(1, 2, 3);
	assert(wc.dot(wc) == 14 && wc != wc + wc);

	PL a(1, 2), b(10, 20);
	a += b; assert(a == PL(11, 22));
	(a -= b) -= b; assert(a == PL(-9, -18));
	assert(a != b && !(a != a));
	vector<PL> v = {PL(0, 0), PL(5, 7)};
	for (auto& t : v) t += b; // as used by MinkowskiSum.h
	assert(v[0] == b && v[1] == PL(15, 27));
	QL p(1, 2, 3), q(1, 2, 4);
	assert(p != q && !(p != p) && p == p);
	assert(mulDiv(4000000000LL, 8000000000000000000LL,
		4000000000000000000LL) == 8000000000LL);
	assert(mulDiv(-7LL, 3LL, 2LL) == -10 && mulDiv(7, 3, 2) == 10);
	assert(mulDiv(1000000000, 2000000000, 1000000000) == 2000000000);
	assert(mulDiv(3.0, 0.5, 4.0) == 0.375);
}

void integers() {
	PD u = PL(3, 4).unit();
	assert(u.x == 0.6 && u.y == 0.8); // used to be (0,0)
	PD n = PL(3, 4).normal();
	assert(n.x == -0.8 && n.y == 0.6);
	PD w = PL(1, 1).unit();
	assert(abs(w.x - sqrt(.5)) < 1e-15 && abs(w.y - sqrt(.5)) < 1e-15);
	PD r = PL(1000, 0).rotate(0.5); // used to truncate to (0,0)
	assert(abs(r.x - 1000*cos(0.5)) < 1e-12);
	assert(abs(r.y - 1000*sin(0.5)) < 1e-12);
	assert(PL(3, 4).dist() == 5 && PL(0, 0).dist() == 0);
	assert(abs(PL(-1, 0).angle() - acos(-1)) < 1e-15);
	// top of the range: dist2 = 8e18 is still exact in long long
	ll big = 2000000000;
	assert(PL(big, -big).dist2() == 8000000000000000000LL);
	assert(PL(big, -big).dist() == sqrt(8e18));
	assert(PL(3*big/5, -4*big/5).dist() == 2e9);
	QD q = QL(3, 4, 12).unit(); // used to be (0,0,0)
	assert(abs(q.x - 3/13.) < 1e-15 && abs(q.y - 4/13.) < 1e-15);
	assert(abs(q.z - 12/13.) < 1e-15);
	assert(QL(3, 4, 12).dist() == 13);
	QD m = QL(2, 0, 0).normal(QL(0, 3, 0));
	assert(m.x == 0 && m.y == 0 && m.z == 1);
	QD t = QL(1000, 0, 0).rotate(0.5, QL(0, 0, 7));
	assert(abs(t.x - 1000*cos(0.5)) < 1e-12);
	assert(abs(t.y - 1000*sin(0.5)) < 1e-12 && t.z == 0);
	big = 1500000000;
	assert(QL(big, -big, big).dist() == sqrt(6.75e18));
	assert(abs(QL(0, 0, -5).theta() - acos(-1)) < 1e-15);
}

void extended() {
	const ld eps = numeric_limits<ld>::epsilon();
	assert(abs(PE(1, 1).dist() - sqrtl(2)) <= 2*eps);
	PE u = PE(1, 1).unit();
	assert(abs(u.x*u.x + u.y*u.y - 1) <= 4*eps);
	assert(abs(PE(1, 1).angle() - atanl(1)) <= eps);
	assert(abs(PE(1e-200L, 0).dist() / 1e-200L - 1) <= 4*eps);
	assert(abs(PE(1e200L, 0).dist() / 1e200L - 1) <= 4*eps);
	assert(abs(PE(3e3000L, 4e3000L).dist() / 5e3000L - 1) <= 4*eps);
	assert(abs(PE(3e-3000L, 4e-3000L).dist() / 5e-3000L - 1) <= 4*eps);
	ld a = 1 / 3.0L;
	PE r = PE(1, 0).rotate(a);
	assert(abs(r.x - cosl(a)) <= eps && abs(r.y - sinl(a)) <= eps);
	assert(abs(QE(1, 1, 1).dist() - sqrtl(3)) <= 2*eps);
	QE q = QE(1, 2, 2).unit();
	assert(abs(q.x - 1/3.0L) <= eps && abs(q.z - 2/3.0L) <= eps);
	assert(abs(QE(1e-2000L, 0, 0).dist() / 1e-2000L - 1) <= 4*eps);
	assert(abs(QE(0, 3e2000L, 4e2000L).dist() / 5e2000L - 1) <= 4*eps);
	QE s = QE(1, 0, 0).rotate(a, QE(0, 0, 2));
	assert(abs(s.x - cosl(a)) <= 2*eps && abs(s.y - sinl(a)) <= 2*eps);
	assert(PE(0, 0).dist() == 0 && QE(0, 0, 0).dist() == 0);
	assert(PD(0, 0).dist() == 0 && QD(0, 0, 0).dist() == 0);
	assert(PD(0, -3e-170).dist() == 3e-170 && PD(5e-324, 0).dist() == 5e-324);
	// double: no overflow or underflow either
	assert(abs(PD(3e200, 4e200).dist() / 5e200 - 1) < 1e-15);
	assert(abs(PD(3e-200, 4e-200).dist() / 5e-200 - 1) < 1e-15);
	PD v = PD(3e200, 4e200).unit();
	assert(abs(v.x - 0.6) < 1e-15 && abs(v.y - 0.8) < 1e-15);
	assert(abs(QD(1e200, 2e200, 2e200).dist() / 3e200 - 1) < 1e-15);
	assert(abs(QD(1e-200, 2e-200, 2e-200).dist() / 3e-200 - 1) < 1e-15);
}

double rnd(mt19937_64& rng, double m) {
	return uniform_real_distribution<double>(-m, m)(rng); }

// double behaviour against long double references
void randomDouble() {
	mt19937_64 rng(12345);
	for (double m : {1e-3, 1., 1e3, 1e9}) rep(it,0,200000) {
		PD p(rnd(rng, m), rnd(rng, m));
		PE e(p.x, p.y);
		ld d = sqrtl(e.x*e.x + e.y*e.y);
		double tol = 4e-16;
		assert(abs(p.dist() - d) <= tol * d);
		PD u = p.unit();
		assert(abs(u.x - e.x/d) <= tol && abs(u.y - e.y/d) <= tol);
		PD n = p.normal();
		assert(abs(n.x + e.y/d) <= tol && abs(n.y - e.x/d) <= tol);
		double a = rnd(rng, 7);
		PD r = p.rotate(a);
		ld c = cosl(a), s = sinl(a);
		assert(abs(r.x - (e.x*c - e.y*s)) <= tol * d);
		assert(abs(r.y - (e.x*s + e.y*c)) <= tol * d);
		assert(p.angle() == atan2(p.y, p.x));
		PL q((ll)rnd(rng, 1e9), (ll)rnd(rng, 1e9));
		if (q == PL(0, 0)) continue;
		assert(abs(q.dist() - sqrtl((ld)q.dist2())) <=
			tol * q.dist());
		PD qu = q.unit();
		assert(abs(qu.dist() - 1) <= tol);
		assert(abs(qu.cross(PD((double)q.x, (double)q.y))) <=
			tol * q.dist());

		QD x(rnd(rng, m), rnd(rng, m), rnd(rng, m));
		QD ax(rnd(rng, m), rnd(rng, m), rnd(rng, m));
		ld xd = sqrtl((ld)x.x*x.x + (ld)x.y*x.y + (ld)x.z*x.z);
		assert(abs(x.dist() - xd) <= tol * xd);
		assert(abs(x.unit().dist() - 1) <= 2 * tol);
		assert(abs(x.theta() - acosl(x.z / xd)) <= 1e-8);
		// Rodrigues with the opposite sign convention check:
		// rotating by a and then by -a is the identity, the
		// length and the axis component are preserved.
		QD y = x.rotate(a, ax), z = y.rotate(-a, ax);
		assert((z - x).dist() <= 8 * tol * xd);
		assert(abs(y.dist() - xd) <= 8 * tol * xd);
		assert(abs(y.dot(ax.unit()) - x.dot(ax.unit())) <= 8*tol*xd);
	}
	// ccw about +z maps +x to +y
	QD y = QD(1, 0, 0).rotate(acos(-1) / 2, QD(0, 0, 1));
	assert(abs(y.x) < 1e-15 && abs(y.y - 1) < 1e-15 && y.z == 0);
}

int main() {
	types();
	operators();
	integers();
	extended();
	randomDouble();
	cout << "Tests passed!" << endl;
}
