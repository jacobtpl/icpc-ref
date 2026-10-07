#include "../utilities/template.h"

#include "../../content/geometry/linearTransformation.h"

typedef complex<long double> C;
mt19937_64 rng(31337);

bool close(P a, C b, long double tol) {
	return fabsl(a.x - b.real()) <= tol && fabsl(a.y - b.imag()) <= tol;
}

int main() {
	// integer grid: oracle is complex arithmetic z -> q0 + (z-p0)*(q1-q0)/(p1-p0)
	const int G = 2;
	rep(ax,-G,G+1) rep(ay,-G,G+1) rep(bx,-G,G+1) rep(by,-G,G+1)
	rep(cx,-G,G+1) rep(cy,-G,G+1) rep(dx,-G,G+1) rep(dy,-G,G+1) {
		if (ax == bx && ay == by) continue;
		P p0(ax, ay), p1(bx, by), q0(cx, cy), q1(dx, dy);
		C a(ax, ay), b(bx, by), c(cx, cy), d(dx, dy);
		assert(close(linearTransformation(p0, p1, q0, q1, p0), c, 1e-12));
		assert(close(linearTransformation(p0, p1, q0, q1, p1), d, 1e-12));
		rep(rx,-3,4) rep(ry,-3,4) {
			C z(rx, ry);
			assert(close(linearTransformation(p0, p1, q0, q1, P(rx, ry)),
				c + (z - a) * (d - c) / (b - a), 1e-11));
		}
	}
	// random doubles at several magnitudes
	for (double r : {1e-3, 1.0, 1e3, 1e6}) {
		uniform_real_distribution<double> U(-r, r);
		rep(it,0,500000) {
			P p0(U(rng), U(rng)), p1(U(rng), U(rng)), q0(U(rng), U(rng)), q1(U(rng), U(rng));
			if ((p1 - p0).dist() < 1e-3 * r) continue; // ill-conditioned
			C a(p0.x, p0.y), b(p1.x, p1.y), c(q0.x, q0.y), d(q1.x, q1.y);
			P r0(U(rng), U(rng)), r1(U(rng), U(rng));
			C z0(r0.x, r0.y), z1(r1.x, r1.y);
			long double tol = 1e-9 * r;
			P t0 = linearTransformation(p0, p1, q0, q1, r0);
			P t1 = linearTransformation(p0, p1, q0, q1, r1);
			assert(close(t0, c + (z0 - a) * (d - c) / (b - a), tol));
			assert(close(t1, c + (z1 - a) * (d - c) / (b - a), tol));
			assert(close(linearTransformation(p0, p1, q0, q1, p0), c, tol));
			assert(close(linearTransformation(p0, p1, q0, q1, p1), d, tol));
			// similarity: all distances scale by |q1-q0|/|p1-p0|
			double k = (q1 - q0).dist() / (p1 - p0).dist();
			assert(abs((t1 - t0).dist() - k * (r1 - r0).dist()) <= 1e-6 * r);
			// inverse map takes us back
			if ((q1 - q0).dist() < 1e-3 * r) continue;
			assert(close(linearTransformation(q0, q1, p0, p1, t0), z0, 1e-6 * r));
		}
	}
	// identity, pure translation, rotation by 90 degrees, scaling
	assert(close(linearTransformation(P(1,2), P(5,-3), P(1,2), P(5,-3), P(7,7)), C(7,7), 1e-12));
	assert(close(linearTransformation(P(0,0), P(1,0), P(10,20), P(11,20), P(3,4)), C(13,24), 1e-12));
	assert(close(linearTransformation(P(0,0), P(1,0), P(0,0), P(0,1), P(3,4)), C(-4,3), 1e-12));
	assert(close(linearTransformation(P(0,0), P(1,0), P(0,0), P(2,0), P(3,4)), C(6,8), 1e-12));
	cout<<"Tests passed!"<<endl;
}
