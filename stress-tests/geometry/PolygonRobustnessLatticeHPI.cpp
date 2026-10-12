#include "../utilities/template.h"
#include "../../content/geometry/HalfPlaneIntersection.h"

int main() {
	using P = Point<double>;
	using IP = Point<ll>;
	mt19937_64 rng(89173);
	vector<ll> sizes = {10, 1000, 100000, 10000000, 1000000000};
	rep(it,0,2000) sizes.push_back(10 + (ll)(rng() % 999999990));
	for (ll n : sizes) {
		vector<IP> tri = {IP(), IP(n-1, n-2), IP(n, n-1)};
		assert(tri[0].cross(tri[1], tri[2]) == 1);
		for (int reflection : {-1, 1}) {
			vector<Ray> rays;
			vector<P> want;
			rep(i,0,3) {
				IP p = tri[i] * reflection, d = (tri[(i+1)%3] - tri[i]) * reflection;
				rays.push_back({P((double)p.x, (double)p.y), P((double)d.x, (double)d.y)});
				want.emplace_back((double)p.x, (double)p.y);
			}
			shuffle(all(rays), rng);
			auto got = halfPlaneIsect(rays);
			sort(all(got)); sort(all(want));
			assert(got == want);
			for (auto& r : rays) r.dp = r.dp * -1;
			assert(halfPlaneIsect(rays).empty());
		}
	}
	cout << "Tests passed!" << endl;
}
