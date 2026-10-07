#include "../utilities/template.h"

#include "../../content/geometry/CircleIntersection.h"

int main() {
	cin.sync_with_stdio(0); cin.tie(0);
	cin.exceptions(cin.failbit);
	srand(2);
	rep(it,0,100000) {
		double rnd[6];
		rep(i,0,6)
			rnd[i] = rand() % 21 - 10;
		P a(rnd[0], rnd[1]);
		P b(rnd[2], rnd[3]);
		double ra = rand() % 10;
		double rb = rand() % 10;
		if (a == b) continue;
		pair<P, P> out;
		bool ret = circleInter(a, b, ra, rb, &out);
		if (ret) {
			assert(abs((out.first - a).dist() - ra) < 1e-9);
			assert(abs((out.second - a).dist() - ra) < 1e-9);
			assert(abs((out.first - b).dist() - rb) < 1e-9);
			assert(abs((out.second - b).dist() - rb) < 1e-9);
		}

		// Hill-climb the answer
		auto func = [&](P x) {
			double d1 = (x - a).dist() - ra;
			double d2 = (x - b).dist() - rb;
			return d1*d1 + d2*d2;
		};
		P start = (a + b) / 2 + (a - b).perp();
		pair<double, P> cur(func(start), start);
		for (double jmp = 100; jmp > 1e-20; jmp /= 2) {
			int iters = 0;
			for (int imp = 1; imp--;) {
				if (++iters == 100) goto skip;
				rep(dx,-1,2) rep(dy,-1,2) {
					P p = cur.second;
					p.x += dx*jmp;
					p.y += dy*jmp;
					pair<double, P> np{func(p), p};
					if (np < cur) cur = np, imp = 1;
				}
			}
		}

		if (abs((cur.second - a).dist() - ra) < 1e-9 &&
		    abs((cur.second - b).dist() - rb) < 1e-9) {
			assert(ret);
			assert((out.first - cur.second).dist() < 1e-6 || (out.second - cur.second).dist() < 1e-6);
		} else {
			assert(!ret);
		}

		// cerr << '.';
		continue;
skip:;
		// Sometimes hill-climbing is slow, for some reason. :(
		// cerr << '#';
	}
	// Constructed instances over several magnitudes: pick the intersection
	// point x first, so the answer is known.
	mt19937_64 gen(5);
	auto uni = [&](double lo, double hi) {
		return lo + (hi - lo) * ((double)(gen() >> 11) / 9007199254740992.0);
	};
	for (double sc : {1e-3, 1.0, 1e3, 1e6}) rep(it,0,200000) {
		P a(uni(-sc, sc), uni(-sc, sc)), b(uni(-sc, sc), uni(-sc, sc));
		P x(uni(-sc, sc), uni(-sc, sc));
		double d = (b - a).dist(), h = (b - a).cross(x - a) / d;
		if (d < 1e-3 * sc || abs(h) < 1e-3 * sc) continue; // ill-conditioned
		double ra = (x - a).dist(), rb = (x - b).dist();
		pair<P, P> out;
		assert(circleInter(a, b, ra, rb, &out));
		// first is to the left of a->b, second is its mirror image
		P y = x - (b - a).perp() * (2 * h / d);
		if (h < 0) swap(x, y);
		assert((out.first - x).dist() < 1e-9 * sc);
		assert((out.second - y).dist() < 1e-9 * sc);
		// clearly disjoint / clearly nested circles
		double r1 = uni(0, d), r2 = uni(0, d - r1) * 0.999;
		assert(!circleInter(a, b, r1, r2, &out));
		assert(!circleInter(a, b, r1 + d * 1.001 + r2, r2, &out));
		assert(!circleInter(a, b, r2, r1 + d * 1.001 + r2, &out));
	}
	// Exact tangency (integer 3-4-5 configurations) and zero radii.
	// (tangency is sqrt-conditioned: error ~ sqrt(eps) * radius, radii <= 1e4)
	auto mx = [&](pair<P, P>& o, P t) {
		return max((o.first - t).dist(), (o.second - t).dist()); };
	const double TOL = 1e-3;
	rep(it,0,200000) {
		int k = (int)(gen() % 1000) + 1, s = (int)(gen() % (5 * k + 1));
		P a((double)((int)(gen() % 2001) - 1000), (double)((int)(gen() % 2001) - 1000));
		P dir((gen() & 1 ? 3 : -3) * k, (gen() & 1 ? 4 : -4) * k);
		if (gen() & 1) swap(dir.x, dir.y);
		P b = a + dir;
		pair<P, P> out;
		// external: r1 + r2 = d, touching at a + dir * s / (5k)
		assert(circleInter(a, b, s, 5 * k - s, &out));
		P t = a + dir * (s / (5.0 * k));
		assert(mx(out, t) < TOL);
		// internal: r1 - r2 = d, touching at a + dir * r1 / d
		assert(circleInter(a, b, 5 * k + s, s, &out));
		t = a + dir * ((5.0 * k + s) / (5.0 * k));
		assert(mx(out, t) < TOL);
		assert(circleInter(b, a, s, 5 * k + s, &out));
		assert(mx(out, t) < TOL);
		// concentric circles with different radii never intersect
		assert(!circleInter(a, a, s, s + 1, &out));
	}
	cout<<"Tests passed!"<<endl;
}
