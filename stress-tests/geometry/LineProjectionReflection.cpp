#include "../utilities/template.h"

#include "../../content/geometry/LineProjRefl.h"
#include "../../content/geometry/lineDistance.h"

typedef Point<double> P;
typedef Point<ll> PL;
int main() {
	cin.sync_with_stdio(0);
	cin.tie(0);
	const int lim = 5;
	for (int i = 0; i < 100000; i++) {
		P p = P(rand() % lim, rand() % lim);
		P a = P(rand() % lim, rand() % lim);
		P b = P(rand() % lim, rand() % lim);
		while (a == b)
			b = P(rand() % lim, rand() % lim);
		auto proj = lineProj(a, b, p, false);
		auto refl = lineProj(a, b, p, true);
		assert(lineDist(a, b, proj) < 1e-8);
		auto manProj = (refl + p) / 2;
		assert((proj-manProj).dist() < 1e-8);
	}
	// Independent checks, with negative and large coordinates: the projection
	// is on the line and p-proj is orthogonal to it; the reflection is at the
	// same distance from a and b as p, on the other side.
	for (int i = 0; i < 500000; i++) {
		int C = i % 3 == 0 ? 3 : i % 3 == 1 ? 1000 : 1000000;
		auto r = [&]() { return (double)(rand() % (2*C+1) - C); };
		P p(r(), r()), a(r(), r()), b(r(), r());
		if (a == b) continue;
		P v = b - a;
		double len = v.dist(), tol = 1e-9 * (C + len);
		auto proj = lineProj(a, b, p), refl = lineProj(a, b, p, true);
		assert(abs(v.cross(proj - a)) / len < tol);
		assert(abs(v.dot(p - proj)) / len < tol);
		assert(abs(lineDist(a, b, proj)) < tol);
		assert(abs(lineDist(a, b, refl) + lineDist(a, b, p)) < tol);
		assert(abs(v.dot(refl - p)) / len < tol);
		assert(abs((refl - a).dist() - (p - a).dist()) < tol);
		assert(abs((refl - b).dist() - (p - b).dist()) < tol);
		// brute force: no point on the line is closer to p than proj
		double best = (p - proj).dist();
		for (int k = -20; k <= 20; k++) {
			P q = proj + v * (k / 7.0);
			assert((p - q).dist() > best - tol);
		}
		// reversing the line gives the same answer
		assert((lineProj(b, a, p) - proj).dist() < tol);
		assert((lineProj(b, a, p, true) - refl).dist() < tol);
	}
	// Integer points: exact whenever the answer is an integer point
	// (axis-parallel lines, and reflections across diagonals).
	for (int i = 0; i < 500000; i++) {
		ll C = i % 2 ? 10 : 1000000;
		auto r = [&]() { return (ll)(rand() % (2*C+1) - C); };
		PL p(r(), r()), a(r(), r());
		ll k = r(); if (!k) k = 1;
		assert(lineProj(a, a + PL(k, 0), p) == PL(p.x, a.y));
		assert(lineProj(a, a + PL(0, k), p) == PL(a.x, p.y));
		assert(lineProj(a, a + PL(k, 0), p, true) == PL(p.x, 2*a.y - p.y));
		assert(lineProj(a, a + PL(0, k), p, true) == PL(2*a.x - p.x, p.y));
		PL d = p - a;
		assert(lineProj(a, a + PL(k, k), p, true) == a + PL(d.y, d.x));
		assert(lineProj(a, a + PL(k, -k), p, true) == a - PL(d.y, d.x));
	}
	cout<<"Tests passed!"<<endl;
}
