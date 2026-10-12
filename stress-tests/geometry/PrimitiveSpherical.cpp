#include "../utilities/template.h"

#include "../../content/geometry/sphericalDistance.h"

typedef long double ld;

// angle between unit vectors as atan2(|u x v|, u . v), long double
ld ref(ld f1, ld t1, ld f2, ld t2) {
	ld a[3] = {sinl(t1)*cosl(f1), sinl(t1)*sinl(f1), cosl(t1)};
	ld b[3] = {sinl(t2)*cosl(f2), sinl(t2)*sinl(f2), cosl(t2)};
	ld c[3] = {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0]};
	return atan2l(sqrtl(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]),
		a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

int main() {
	const double PI = acos(-1);
	const ld PIL = acosl(-1);
	// reviewed reproducer: asin(d/2) with d > 2 gave NaN
	double r = sphericalDistance(0.008, 0.96, 0.008 + PI, PI - 0.96, 1);
	assert(abs(r - PI) < 1e-15);
	// exactly antipodal grid
	rep(i,1,2000) rep(j,1,100) {
		double f = i * 0.001, t = j * 0.03;
		double d = sphericalDistance(f, t, f + PI, PI - t, 1.0);
		assert(!std::isnan(d) && abs(d - PI) < 1e-15);
		assert(d <= PI);
	}
	// near antipodal: used to lose half of the digits
	for (double delta : {1e-3, 1e-6, 1e-8, 1e-10, 1e-13}) {
		double d = sphericalDistance(0, PI / 2, PI - delta, PI / 2, 1);
		assert(abs(d - (PIL - delta)) < 1e-15);
		d = sphericalDistance(1, 1, 1 + PI, PI - 1 + delta, 1);
		assert(abs(d - (PIL - delta)) < 1e-15);
	}
	// nearby points keep their relative accuracy on the equator
	for (double delta : {1e-3, 1e-6, 1e-9, 1e-12})
		assert(abs(sphericalDistance(0, PI / 2, delta, PI / 2, 1) / delta
			- 1) < 1e-15);
	assert(sphericalDistance(0.3, 0.7, 0.3, 0.7, 5) == 0);
	assert(abs(sphericalDistance(0, 0, 0, PI / 2, 2) - PI) < 1e-15);
	assert(abs(sphericalDistance(0, PI / 2, PI / 2, PI / 2, 3) - 1.5 * PI)
		< 1e-15);

	mt19937_64 rng(31337);
	uniform_real_distribution<double> lon(-PI, PI), lat(0, PI);
	double worst = 0;
	rep(it,0,1000000) {
		double f1 = lon(rng), t1 = lat(rng), f2 = lon(rng), t2 = lat(rng);
		if (it % 4 == 1) f2 = f1 + PI, t2 = PI - t1; // antipodal
		if (it % 4 == 2) f2 = f1 + PI + lon(rng) * 1e-7, t2 = PI - t1;
		if (it % 4 == 3) f2 = f1 + lon(rng) * 1e-7, t2 = t1 + 1e-9;
		double d = sphericalDistance(f1, t1, f2, t2, 1);
		assert(!std::isnan(d) && d >= 0 && d <= PI);
		assert(d == sphericalDistance(f2, t2, f1, t1, 1));
		worst = max(worst, (double)abs(d - ref(f1, t1, f2, t2)));
		double R = 6371;
		assert(abs(sphericalDistance(f1, t1, f2, t2, R) - R * d)
			<= 1e-15 * R * PI);
	}
	assert(worst < 1e-15);
	cout << "Tests passed!" << endl;
}
