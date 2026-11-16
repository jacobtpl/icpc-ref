#include "../utilities/template.h"
#include "../../content/geometry/SegmentIntersection.h"

using IP = Point<ll>;
using DP = Point<double>;
using LP = Point<long double>;
using FP = Point<Fraction>;
typedef __int128 LL;

int main() {
	mt19937_64 rng(832704);
	for (ll scale : {1LL, 1000000LL, 1000000000LL}) {
		IP a(-scale, -scale), b(scale, scale);
		IP c(-scale, scale), d(scale, -scale);
		assert(segInter(a, b, c, d) == vector<LP>{LP(0, 0)});
		assert(segInter<Fraction>(a, b, c, d) == vector<FP>{FP(0, 0)});
	}
	rep(it,0,100000) {
		ll R = it & 1 ? 1000000000 : 1000000;
		auto r = [&]() { return (ll)(rng() % (2 * R + 1)) - R; };
		IP a(r(), r()), b(r(), r()), c(r(), r()), d(r(), r());
		LL oa = c.cross(d, a), ob = c.cross(d, b),
			oc = a.cross(b, c), od = a.cross(b, d);
		if (sgn(oa) * sgn(ob) >= 0 || sgn(oc) * sgn(od) >= 0) continue;
		LL nx = a.x * ob - b.x * oa, ny = a.y * ob - b.y * oa, den = ob - oa;
		for (auto res : {segInter(a, b, c, d), segInter(b, a, d, c)}) {
			assert(sz(res) == 1);
			assert(fabsl(res[0].x - (long double)nx / (long double)den) < 1e-6);
			assert(fabsl(res[0].y - (long double)ny / (long double)den) < 1e-6);
		}
		if (R == 1000000) { // exact rational answer fits in a Fraction
			FP want(Fraction::fromWide(nx, den), Fraction::fromWide(ny, den));
			assert(segInter<Fraction>(a, b, c, d) == vector<FP>{want});
			assert(segInter<Fraction>(b, a, d, c) == vector<FP>{want});
		}
	}
	for (double offset : {0., 1e9, 1e12}) {
		DP a(offset, offset), b(offset + 2, offset + 2);
		DP c(offset, offset + 2), d(offset + 2, offset);
		assert(segInter(a, b, c, d) == vector<DP>{DP(offset+1, offset+1)});
		assert(segInter(b, a, d, c) == vector<DP>{DP(offset+1, offset+1)});
	}
	cout << "Tests passed!" << endl;
}
