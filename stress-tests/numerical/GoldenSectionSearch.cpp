#include "../utilities/template.h"

#include "../../content/numerical/GoldenSectionSearch.h"

// gss takes a plain function pointer, so the test functions are
// parameterised through globals.
double C, K; int evals;
double fQuad(double x) { evals++; return K * (x - C) * (x - C); }
double fAbs(double x) { evals++; return K * fabs(x - C); }
double fQuart(double x) { evals++; double d = x - C; return d * d * d * d; }
double fAsym(double x) { evals++; double d = x - C; return d < 0 ? -K * d : 5 * d * d + d; }
double fCosh(double x) { evals++; return cosh((x - C) / 50); }
double fInc(double x) { evals++; return K * x; }
double fDec(double x) { evals++; return -K * x; }
double fEx(double x) { return 4+x+.3*x*x; }

mt19937_64 rng(12345);
double rnd(double lo, double hi) {
	return uniform_real_distribution<double>(lo, hi)(rng);
}

void check(double a, double b, double (*f)(double), double want, double tol) {
	evals = 0;
	double x = gss(a, b, f);
	if (!(a <= x && x <= b) || fabs(x - want) > tol) {
		fprintf(stderr, "gss(%.17g, %.17g) = %.17g, want %.17g\n", a, b, x, want);
		abort();
	}
	// O(log((b-a)/eps)) evaluations, base 1/r
	double bound = 3 + max(0.0, log((b - a) / 1e-7) / log(1.6180339));
	assert(evals <= bound + 2);
}

int main() {
	// header usage example: 4+x+.3x^2 has its minimum at -5/3
	assert(fabs(gss(-1000, 1000, fEx) + 5.0/3) < 1e-6);

	double (*fs[])(double) = {fQuad, fAbs, fQuart, fAsym, fCosh};
	rep(it,0,300000) {
		double scale = pow(10, rnd(-3, 4));
		double a = rnd(-scale, scale), b = rnd(-scale, scale);
		if (a > b) swap(a, b);
		C = rnd(a, b); K = pow(10, rnd(-2, 2));
		int which = (int)(rng() % 5);
		// quartic and cosh are numerically flat around the minimum (ties
		// in f), so only a loose tolerance is meaningful for them
		double tol = which == 2 ? 2e-3 : which == 4 ? 1e-4 : 2e-7;
		check(a, b, fs[which], C, tol);
		// minimum outside of the interval: answer is an endpoint
		C = b + rnd(0, scale); check(a, b, fs[which], b, tol);
		C = a - rnd(0, scale); check(a, b, fs[which], a, tol);
		check(a, b, fInc, a, 2e-7);
		check(a, b, fDec, b, 2e-7);
		// minimum exactly on an endpoint
		C = a; check(a, b, fs[which % 2], a, 2e-7);
		C = b; check(a, b, fs[which % 2], b, 2e-7);
	}
	// degenerate intervals
	rep(it,0,1000) {
		double a = rnd(-1e6, 1e6);
		C = a; K = 1;
		check(a, a, fQuad, a, 0);
		check(a, a + 1e-9, fQuad, a, 1e-9);
		check(a, a + 5e-8, fAbs, a, 5e-8);
	}
	// Largest magnitudes that are safe: the spacing of doubles must stay
	// below eps = 1e-7, i.e. |x| < 2^29 ~ 5.3e8.
	rep(it,0,20000) {
		double a = rnd(-5e8, 5e8), b = rnd(-5e8, 5e8);
		if (a > b) swap(a, b);
		C = rnd(a, b); K = 1;
		check(a, b, it % 2 ? fQuad : fAbs, C, 3e-7);
	}
#ifdef GSS_LARGE
	// Opt-in: for |x| >= 2^29 adjacent doubles are more than eps = 1e-7
	// apart, so b-a never drops below eps and gss never terminates.
	{
		static int cnt;
		auto f = [](double x) {
			if (++cnt > 1000000) {
				puts("gss(0, 2e9) did not terminate after 1e6 evaluations");
				exit(1);
			}
			return fabs(x - 1e9);
		};
		double x = gss(0, 2e9, f);
		assert(fabs(x - 1e9) < 1e-6);
	}
#endif
	cout<<"Tests passed!"<<endl;
}
