#include "../utilities/template.h"

#include "../../content/geometry/sphericalDistance.h"

typedef long double ld;
mt19937_64 rng(271828);
const double PI = acos(-1);

// oracle: angle between the unit vectors via atan2(|a x b|, a . b)
ld oracle(ld f1, ld t1, ld f2, ld t2, ld radius) {
	ld a[3] = {sinl(t1)*cosl(f1), sinl(t1)*sinl(f1), cosl(t1)};
	ld b[3] = {sinl(t2)*cosl(f2), sinl(t2)*sinl(f2), cosl(t2)};
	ld c[3] = {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
	ld dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
	return radius * atan2l(sqrtl(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]), dot);
}

void test(double f1, double t1, double f2, double t2, double r) {
	double got = sphericalDistance(f1, t1, f2, t2, r);
	ld want = oracle(f1, t1, f2, t2, r);
	assert(!std::isnan(got));
	// asin(d/2) is ill-conditioned near antipodal points (error ~ sqrt(eps))
	assert(fabsl(got - want) <= 1e-7 * r);
	assert(got >= 0 && got <= PI * r * (1 + 1e-15));
	assert(sphericalDistance(f2, t2, f1, t1, r) == got);
}

int main() {
	uniform_real_distribution<double> F(-2 * PI, 2 * PI), T(0, PI), R(0.1, 1e4);
	// grid of special angles (poles, equator, wrap-around)
	vector<double> fs, ts;
	rep(i,-8,9) fs.push_back(i * PI / 4);
	rep(i,0,9) ts.push_back(i * PI / 8);
	for (double f1 : fs) for (double t1 : ts) for (double f2 : fs) for (double t2 : ts)
		test(f1, t1, f2, t2, 1), test(f1, t1, f2, t2, 6371);
	// random points
	rep(it,0,1000000) test(F(rng), T(rng), F(rng), T(rng), R(rng));
	// nearly identical points
	rep(it,0,300000) {
		double f = F(rng), t = T(rng), e = pow(10, -uniform_real_distribution<double>(1, 12)(rng));
		double f2 = f + e * F(rng), t2 = min(PI, max(0.0, t + e * F(rng))), r = R(rng);
		double got = sphericalDistance(f, t, f2, t2, r);
		ld want = oracle(f, t, f2, t2, r);
		assert(fabsl(got - want) <= 1e-12 * r);
	}
	// same point
	rep(it,0,100000) {
		double f = F(rng), t = T(rng);
		assert(sphericalDistance(f, t, f, t, R(rng)) == 0);
	}
	// antipodal points: distance must be pi*radius, never nan
	rep(it,0,1000000) {
		double f = F(rng), t = T(rng), r = R(rng);
		test(f, t, f + PI, PI - t, r);
		test(f, t, f - PI, PI - t, r);
	}
	// known values
	assert(abs(sphericalDistance(0, 0, 0, PI, 1) - PI) < 1e-9);
	assert(abs(sphericalDistance(0, PI/2, PI/2, PI/2, 2) - PI) < 1e-9);
	assert(abs(sphericalDistance(1, 0, 2, PI/2, 3) - 3*PI/2) < 1e-9);
	cout<<"Tests passed!"<<endl;
}
