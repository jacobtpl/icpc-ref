#include "../utilities/template.h"

#include "../../content/numerical/Integrate.h"

mt19937_64 rng(777);
double rnd(double lo, double hi) {
	return uniform_real_distribution<double>(lo, hi)(rng);
}

void close(double got, double want, double rel, const char* what) {
	if (!(fabs(got - want) <= rel * max(1.0, fabs(want)))) {
		fprintf(stderr, "%s: got %.15g want %.15g\n", what, got, want);
		abort();
	}
}

int main() {
	// Simpson's rule is exact for cubics, for every n >= 1
	rep(it,0,100000) {
		double c[4];
		rep(i,0,4) c[i] = rnd(-5, 5);
		double a = rnd(-10, 10), b = rnd(-10, 10);
		if (it % 50 == 0) b = a;
		auto f = [&](double x) { return ((c[3]*x + c[2])*x + c[1])*x + c[0]; };
		auto F = [&](double x) {
			return (((c[3]/4*x + c[2]/3)*x + c[1]/2)*x + c[0])*x; };
		int n = it % 3 == 0 ? 1000 : 1 + (int)(rng() % 40);
		double got = it % 3 == 0 ? quad(a, b, f) : quad(a, b, f, n);
		close(got, F(b) - F(a), 1e-9, "cubic");
		// reversed limits flip the sign
		close(quad(b, a, f, n), F(a) - F(b), 1e-9, "cubic reversed");
	}
	// smooth functions with the default n = 1000: the error is
	// about (b-a) h^4 max|f''''| / 180 with h = (b-a)/2000
	rep(it,0,20000) {
		double a = rnd(-10, 10), b = rnd(-10, 10), w = rnd(0.1, 3);
		close(quad(a, b, [&](double x) { return sin(w * x); }),
			(cos(w * a) - cos(w * b)) / w, 3e-7, "sin");
		close(quad(a, b, [&](double x) { return exp(x / 3); }),
			3 * (exp(b / 3) - exp(a / 3)), 3e-7, "exp");
		close(quad(a, b, [&](double x) { return 1 / (1 + x * x); }),
			atan(b) - atan(a), 3e-7, "atan");
	}
	// error is O(h^4): doubling n divides the error by ~16
	{
		auto f = [](double x) { return exp(x); };
		double ex = exp(1.0) - 1, e1 = fabs(quad(0, 1, f, 4) - ex),
			e2 = fabs(quad(0, 1, f, 8) - ex);
		assert(e1 / e2 > 14 && e1 / e2 < 18);
	}
	// number of evaluations is exactly 2n+1
	{
		int cnt = 0;
		quad(0, 1, [&](double x) { cnt++; return x; }, 123);
		assert(cnt == 247);
	}
	cout<<"Tests passed!"<<endl;
}
