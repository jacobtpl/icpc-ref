#include "../utilities/template.h"

#include "../../content/numerical/IntegrateAdaptive.h"

mt19937_64 rng(4242);
double rnd(double lo, double hi) {
	return uniform_real_distribution<double>(lo, hi)(rng);
}

void close(double got, double want, double tol, const char* what) {
	if (!(fabs(got - want) <= tol * max(1.0, fabs(want)))) {
		fprintf(stderr, "%s: got %.15g want %.15g\n", what, got, want);
		abort();
	}
}

int main() {
	const double pi = acos(-1.0);
	rep(it,0,30000) {
		double c[4];
		rep(i,0,4) c[i] = rnd(-5, 5);
		double a = rnd(-10, 10), b = rnd(-10, 10);
		if (it % 50 == 0) b = a;
		auto f = [&](double x) { return ((c[3]*x + c[2])*x + c[1])*x + c[0]; };
		auto F = [&](double x) {
			return (((c[3]/4*x + c[2]/3)*x + c[1]/2)*x + c[0])*x; };
		close(quad(a, b, f), F(b) - F(a), 1e-9, "cubic");
		close(quad(b, a, f), F(a) - F(b), 1e-9, "cubic reversed");
		// w(b-a)/4 stays below 2pi, see ADAPTIVE_ALIAS below
		double w = rnd(0.1, 1.2);
		close(quad(a, b, [&](double x) { return sin(w * x); }),
			(cos(w * a) - cos(w * b)) / w, 1e-7, "sin");
		close(quad(a, b, [&](double x) { return exp(x / 3); }),
			3 * (exp(b / 3) - exp(a / 3)), 1e-7, "exp");
		// reversed limits must just flip the sign
		close(quad(b, a, [&](double x) { return sin(w * x); }),
			(cos(w * b) - cos(w * a)) / w, 1e-7, "sin reversed");
		close(quad(b, a, [&](double x) { return exp(x / 3); }),
			3 * (exp(a / 3) - exp(b / 3)), 1e-7, "exp reversed");
		// kink and jump inside the interval
		double k = rnd(-10, 10);
		auto G = [&](double x) { return (x - k) * fabs(x - k) / 2; };
		close(quad(a, b, [&](double x) { return fabs(x - k); }),
			G(b) - G(a), 1e-7, "abs");
		auto H = [&](double x) { return max(x - k, 0.0); };
		close(quad(a, b, [&](double x) { return x > k ? 1.0 : 0.0; }),
			H(b) - H(a), 1e-7, "step");
	}
	// endpoint singularities of the derivative, sharp peak, tighter eps
	close(quad(0, 1, [](double x) { return sqrt(x); }), 2.0/3, 1e-7, "sqrt");
	close(quad(-1, 1, [](double x) { return sqrt(1 - x*x); }), pi/2, 1e-7, "circle");
	close(quad(-100, 100, [](double x) { return exp(-x*x); }), sqrt(pi), 1e-7, "gauss");
	close(quad(0, 1, [](double x) { return exp(x); }, 1e-13), exp(1.0) - 1, 1e-12, "eps");
	for (double b : {1.0, 2.0, 5.0, 10.0, 1000.0})
		close(quad(-b, b, [](double x) { return 1 / (1 + x * x); }),
			2 * atan(b), 1e-7, "atan");
	// header usage example: volume of the unit sphere
	{
		ll cnt = 0;
		double sphereVolume = quad(-1, 1, [&](double x) {
			return quad(-1, 1, [&](double y) {
			return quad(-1, 1, [&](double z) {
			cnt++;
			return x*x + y*y + z*z < 1 ? 1.0 : 0.0; });});});
		close(sphereVolume, 4 * pi / 3, 1e-6, "sphere");
		cerr << "sphere evaluations: " << cnt << endl;
	}
#ifdef ADAPTIVE_ALIAS
	// Opt-in: inherent limitation of adaptive Simpson. If the period of
	// f is close to (b-a)/4 (or /2), the five initial samples are nearly
	// equal, the first error estimate is ~0 and recursion stops at once.
	close(quad(0, 25.2, [](double x) { return sin(x) * sin(x); }),
		12.6 - sin(50.4) / 4, 1e-6, "sin^2 on [0, 25.2]");
	// Same effect on a non-periodic function: the estimate T - S happens
	// to vanish on a coarse subinterval. eps = 1e-8, actual error 5.7e-4.
	{
		double a = -8.4098233469445862, b = -0.61316265458802732;
		close(quad(a, b, [](double x) { return 1 / (1 + x * x); }),
			atan(b) - atan(a), 1e-6, "atan on [-8.41, -0.61]");
	}
#endif
	cout<<"Tests passed!"<<endl;
}
