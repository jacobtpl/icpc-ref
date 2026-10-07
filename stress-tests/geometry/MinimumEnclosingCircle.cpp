#include "../utilities/template.h"

#include "../../content/geometry/MinimumEnclosingCircle.h"

typedef long double ld;

// brute force: smallest circle through 1, 2 or 3 of the points that encloses
// all of them, in long double. O(n^4).
ld slowMec(const vector<P>& ps) {
	int n = sz(ps);
	ld best = 1e300L;
	auto tryC = [&](ld x, ld y, ld r) {
		if (r >= best) return;
		for (auto &p : ps) if (hypotl(p.x - x, p.y - y) > r + 1e-12L * max((ld)1, r)) return;
		best = r;
	};
	rep(i,0,n) {
		tryC(ps[i].x, ps[i].y, 0);
		rep(j,0,i) {
			ld mx = ((ld)ps[i].x + ps[j].x) / 2, my = ((ld)ps[i].y + ps[j].y) / 2;
			tryC(mx, my, hypotl(ps[i].x - mx, ps[i].y - my));
			rep(k,0,j) {
				ld bx = (ld)ps[j].x - ps[i].x, by = (ld)ps[j].y - ps[i].y;
				ld cx = (ld)ps[k].x - ps[i].x, cy = (ld)ps[k].y - ps[i].y;
				ld d = 2 * (bx * cy - by * cx);
				if (d == 0) continue;
				ld ux = (cy * (bx * bx + by * by) - by * (cx * cx + cy * cy)) / d;
				ld uy = (bx * (cx * cx + cy * cy) - cx * (bx * bx + by * by)) / d;
				tryC(ps[i].x + ux, ps[i].y + uy, hypotl(ux, uy));
			}
		}
	}
	return best;
}

void checkMec(const vector<P>& ps) {
	pair<P, double> pa = mec(ps);
	double r = pa.second;
	ld opt = slowMec(ps);
	ld tol = 1e-7L * max((ld)1e-9, opt);
	if (!(abs(r - opt) <= tol)) {
		cerr << setprecision(17) << "radius " << r << " expected " << (double)opt << " for";
		for (auto &p : ps) cerr << ' ' << p;
		cerr << endl;
		abort();
	}
	for (auto &p : ps) assert((p - pa.first).dist() <= r + (double)tol);
}

void testAgainstBrute() {
	mt19937_64 rng(3);
	auto ri = [&](ll a, ll b) { return (double)uniform_int_distribution<ll>(a, b)(rng); };
	// n = 1, duplicates, collinear
	checkMec({P(3, -4)});
	checkMec({P(3, -4), P(3, -4), P(3, -4)});
	checkMec({P(0, 0), P(1, 1), P(2, 2), P(5, 5), P(-7, -7)});
	checkMec({P(1e9, 1e9), P(-1e9, -1e9), P(1e9, -1e9), P(-1e9, 1e9)});
	// random integer points at many scales, lots of duplicates for small C
	for (ll C : {1LL, 3LL, 10LL, 1000LL, 1000000LL, 1000000000LL}) rep(it,0,30000) {
		int n = (int)ri(1, 9);
		vector<P> ps;
		rep(i,0,n) ps.emplace_back(ri(-C, C), ri(-C, C));
		checkMec(ps);
	}
	// collinear points with large coordinates
	rep(it,0,30000) {
		int n = (int)ri(1, 9);
		double dx = ri(-1000, 1000), dy = ri(-1000, 1000), ox = ri(-1e8, 1e8), oy = ri(-1e8, 1e8);
		vector<P> ps;
		rep(i,0,n) { double t = ri(-100000, 100000); ps.emplace_back(ox + dx * t, oy + dy * t); }
		checkMec(ps);
	}
	// many cocircular lattice points: x^2 + y^2 = 5^2 * 13^2 * 17^2, shifted
	{
		vector<P> circ;
		const int R = 5 * 13 * 17;
		rep(x,-R,R+1) {
			int y = (int)llround(sqrt((double)R * R - (double)x * x));
			if (x * x + y * y == R * R) { circ.emplace_back(x, y); if (y) circ.emplace_back(x, -y); }
		}
		assert(sz(circ) > 30);
		rep(it,0,30000) {
			int n = (int)ri(1, 9);
			double ox = ri(-1e6, 1e6), oy = ri(-1e6, 1e6);
			vector<P> ps;
			rep(i,0,n) ps.push_back(circ[(int)ri(0, sz(circ) - 1)] + P(ox, oy));
			if (ri(0, 1)) ps.emplace_back(ox + ri(-R / 2, R / 2), oy + ri(-R / 2, R / 2));
			checkMec(ps);
		}
	}
	// non-integer points in a tiny cluster far from the origin
	rep(it,0,30000) {
		int n = (int)ri(1, 9);
		double ox = ri(-1000, 1000), oy = ri(-1000, 1000);
		vector<P> ps;
		rep(i,0,n) ps.emplace_back(ox + ri(-1000, 1000) * 1e-3, oy + ri(-1000, 1000) * 1e-3);
		checkMec(ps);
	}
}

int main() {
	testAgainstBrute();
	srand(2);
	rep(it,0,1000000) {
		int N = rand() % 20 + 1;
		// int N = 4;
		vector<P> ps;
		rep(i,0,N) {
			ps.emplace_back(rand() % 21 - 10, rand() % 21 - 10);
		}

		pair<P, double> pa = mec(ps);
		P mid = pa.first;
		double rad = pa.second;
		double maxDist = 0;
		for(auto &p: ps) {
			maxDist = max(maxDist, (p - mid).dist());
		}

		assert(abs(maxDist - rad) < 1e-6);

		rep(it2,0,50) {
			P q2 = mid - P(0, 1e-6).rotate(it2);
			for(auto &p: ps) {
				if((p - q2).dist() > rad - 1e-7) goto fail;
			}
			assert(0);
fail:;
		}
	}
	cout<<"Tests passed!"<<endl;
}
