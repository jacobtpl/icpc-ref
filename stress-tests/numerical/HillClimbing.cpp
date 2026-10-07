#include "../utilities/template.h"

#include "../../content/numerical/HillClimbing.h"

mt19937_64 rng(2024);
double rnd(double lo, double hi) {
	return uniform_real_distribution<double>(lo, hi)(rng);
}

int main() {
	ll evals = 0;
	// convex quadratics (rotated ellipses), minimised at (cx, cy)
	rep(it,0,600) {
		double scale = pow(10, rnd(-2, 6));
		double cx = rnd(-scale, scale), cy = rnd(-scale, scale);
		double th = rnd(0, 7), co = cos(th), si = sin(th);
		// condition number at most 1000, see HILL_ILLCOND below
		double k1 = pow(10, rnd(-1.5, 1.5)), k2 = pow(10, rnd(-1.5, 1.5));
		double off = it % 2 ? rnd(-100, 100) : 0;
		auto f = [&](P p) {
			evals++;
			double x = p[0] - cx, y = p[1] - cy;
			double u = co * x + si * y, v = -si * x + co * y;
			return k1 * u * u + k2 * v * v + off;
		};
		double s2 = pow(10, rnd(-2, 8));
		P start{{rnd(-s2, s2), rnd(-s2, s2)}};
		auto res = hillClimb(start, f);
		// with an offset f only resolves k*d^2 down to ~1e-14, so the
		// position is checked only for off = 0, the value always
		double tol = off != 0 ? 1e100 : 1e-6 * max(1.0, scale);
		if (fabs(res.second[0] - cx) > tol || fabs(res.second[1] - cy) > tol
				|| res.first != f(res.second) || res.first > f(start)
				|| res.first > off + 1e-9) {
			fprintf(stderr, "quad: got (%.12g, %.12g) want (%.12g, %.12g)\n",
				res.second[0], res.second[1], cx, cy);
			fprintf(stderr, "k1=%g k2=%g th=%g off=%g start=(%g, %g) f=%g\n",
				k1, k2, th, off, start[0], start[1], res.first);
			abort();
		}
	}
	// non-smooth unimodal functions: L1, Linf and L2 distance
	rep(it,0,300) {
		double cx = rnd(-1e5, 1e5), cy = rnd(-1e5, 1e5);
		int which = it % 3;
		auto f = [&](P p) {
			double x = fabs(p[0] - cx), y = fabs(p[1] - cy);
			return which == 0 ? x + y : which == 1 ? max(x, y) : hypot(x, y);
		};
		P start{{rnd(-1e7, 1e7), rnd(-1e7, 1e7)}};
		auto res = hillClimb(start, f);
		if (fabs(res.second[0] - cx) > 1e-6 || fabs(res.second[1] - cy) > 1e-6) {
			fprintf(stderr, "norm %d: got (%.12g, %.12g) want (%.12g, %.12g)\n", which,
				res.second[0], res.second[1], cx, cy);
			abort();
		}
	}
	// start exactly at the optimum, constant function, far away start
	{
		auto f = [](P p) { return p[0] * p[0] + p[1] * p[1]; };
		auto res = hillClimb(P{{0, 0}}, f);
		assert(res.first == 0 && res.second == (P{{0, 0}}));
		auto g = [](P) { return 42.0; };
		auto r2 = hillClimb(P{{3, 4}}, g);
		assert(r2.first == 42);
		auto r3 = hillClimb(P{{1e10, -1e10}}, f);
		assert(r3.first < 1e-12);
	}
#ifdef HILL_ILLCOND
	// Opt-in: a convex quadratic whose valley is narrow (condition
	// number 1e4) and neither axis-aligned nor diagonal. The 8 fixed directions stop
	// giving improvements while jmp is still large, and the result is
	// thousands of units away from the minimum at (0, 0).
	{
		auto f = [](P p) {
			double u = 0.8 * p[0] + 0.6 * p[1], v = -0.6 * p[0] + 0.8 * p[1];
			return u * u + 1e4 * v * v;
		};
		auto res = hillClimb(P{{3000, 1000}}, f);
		fprintf(stderr, "ill-conditioned: got (%.6g, %.6g) f=%.6g, want (0, 0)\n",
			res.second[0], res.second[1], res.first);
		assert(fabs(res.second[0]) < 1e-6 && fabs(res.second[1]) < 1e-6);
	}
#endif
	// number of evaluations is fixed: 97 halvings * 100 * 9, plus one
	assert(evals % 600 == 0 && evals / 600 > 87000 && evals / 600 < 88000);
	cout<<"Tests passed!"<<endl;
}
