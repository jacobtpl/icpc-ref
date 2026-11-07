#include "../utilities/template.h"
#include "../../content/geometry/SegmentIntersection.h"

using IP = Point<ll>;
using DP = Point<double>;

int main() {
	mt19937_64 rng(832704);
	for (ll scale : {1LL, 1000000LL, 1000000000LL}) {
		IP a(-scale, -scale), b(scale, scale);
		IP c(-scale, scale), d(scale, -scale);
		assert(segInter(a, b, c, d) == vector<IP>{IP(0, 0)});
	}
	rep(it,0,100000) {
		auto r = [&]() { return (ll)(rng() % 2000000001) - 1000000000; };
		IP a(r(), r()), b(r(), r()), c(r(), r()), d(r(), r());
		__int128 oa = c.cross(d, a), ob = c.cross(d, b),
			oc = a.cross(b, c), od = a.cross(b, d);
		if (sgn(oa) * sgn(ob) >= 0 || sgn(oc) * sgn(od) >= 0) continue;
		IP want((ll)((a.x * ob - b.x * oa) / (ob - oa)),
			(ll)((a.y * ob - b.y * oa) / (ob - oa)));
		assert(segInter(a, b, c, d) == vector<IP>{want});
		assert(segInter(b, a, d, c) == vector<IP>{want});
	}
	for (double offset : {0., 1e9, 1e12}) {
		DP a(offset, offset), b(offset + 2, offset + 2);
		DP c(offset, offset + 2), d(offset + 2, offset);
		assert(segInter(a, b, c, d) == vector<DP>{DP(offset+1, offset+1)});
		assert(segInter(b, a, d, c) == vector<DP>{DP(offset+1, offset+1)});
	}
	cout << "Tests passed!" << endl;
}
