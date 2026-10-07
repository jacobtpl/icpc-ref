#include "../utilities/template.h"

#include "../../content/geometry/Point3D.h"
#include "../../content/geometry/PolyhedronVolume.h"

mt19937_64 rng(13);
ll ri(ll a, ll b) { return uniform_int_distribution<ll>(a, b)(rng); }

struct Tri { int a, b, c; };
typedef Point3D<ll> PL;
typedef Point3D<double> PD;

// unit cube, vertex index = x + 2y + 4z, all faces pointing outwards
const vector<Tri> cubeTris = {
	{0,2,3}, {0,3,1}, {4,5,7}, {4,7,6}, {0,1,5}, {0,5,4},
	{2,6,7}, {2,7,3}, {0,4,6}, {0,6,2}, {1,3,7}, {1,7,5}};

template<class T> vector<Point3D<double>> toDouble(const vector<Point3D<T>>& v) {
	vector<PD> r;
	for (auto& p : v) r.emplace_back((double)p.x, (double)p.y, (double)p.z);
	return r;
}

int main() {
	{
		vector<PD> p;
		rep(i,0,8) p.emplace_back(i & 1, i >> 1 & 1, i >> 2 & 1);
		assert(abs(signedPolyVolume(p, cubeTris) - 1) < 1e-12);
		vector<Tri> flipped = cubeTris;
		for (auto& t : flipped) swap(t.b, t.c); // inward faces: negative
		assert(abs(signedPolyVolume(p, flipped) + 1) < 1e-12);
		assert(signedPolyVolume(p, vector<Tri>()) == 0);
	}
	rep(it,0,300000) {
		ll lim = it % 3 == 0 ? 1000 : it % 3 == 1 ? 3 : 100000;
		PL off(ri(-lim, lim), ri(-lim, lim), ri(-lim, lim));
		{ // affine image of the unit cube: signed volume = det(M)
			ll m[3][3];
			rep(i,0,3) rep(j,0,3) m[i][j] = ri(-lim, lim);
			double det = 0;
			rep(i,0,3) det += (double)m[0][i] * ((double)m[1][(i+1)%3] * (double)m[2][(i+2)%3]
				- (double)m[1][(i+2)%3] * (double)m[2][(i+1)%3]);
			vector<PL> p;
			rep(i,0,8) {
				ll c[3] = {i & 1, i >> 1 & 1, i >> 2 & 1};
				PL q = off;
				rep(j,0,3) q = q + PL(m[0][j], m[1][j], m[2][j]) * c[j];
				p.push_back(q);
			}
			double tol = 1e-9 * (double)lim * (double)lim * (double)lim;
			assert(abs(signedPolyVolume(toDouble(p), cubeTris) - det) <= tol);
		}
		{ // tetrahedron: |det| / 6 after orienting the faces outwards
			vector<PL> p(4);
			for (auto& q : p) q = off + PL(ri(-lim, lim), ri(-lim, lim), ri(-lim, lim));
			PL a = p[1] - p[0], b = p[2] - p[0], c = p[3] - p[0];
			double det = (double)a.x * ((double)b.y*(double)c.z - (double)b.z*(double)c.y)
				- (double)a.y * ((double)b.x*(double)c.z - (double)b.z*(double)c.x)
				+ (double)a.z * ((double)b.x*(double)c.y - (double)b.y*(double)c.x);
			vector<Tri> t = {{0,2,1}, {0,1,3}, {1,2,3}, {0,3,2}};
			if (det < 0) for (auto& f : t) swap(f.b, f.c);
			double tol = 1e-9 * (double)lim * (double)lim * (double)lim;
			assert(abs(signedPolyVolume(toDouble(p), t) - abs(det) / 6) <= tol);
		}
		{ // cone over an arbitrary polygon in the plane z = z0: area * h / 3
			int n = (int)ri(3, 9);
			ll h = ri(1, lim);
			vector<PL> p;
			rep(i,0,n) p.push_back(PL(off.x + ri(-lim, lim), off.y + ri(-lim, lim), off.z));
			double area2 = 0;
			rep(i,0,n) {
				PL u = p[i] - off, v = p[(i+1)%n] - off;
				area2 += (double)u.x * (double)v.y - (double)u.y * (double)v.x;
			}
			p.push_back(PL(off.x + ri(-lim, lim), off.y + ri(-lim, lim), off.z + h));
			vector<Tri> t;
			rep(i,1,n-1) t.push_back({0, i+1, i}); // base, facing down
			rep(i,0,n) t.push_back({i, (i+1)%n, n}); // sides
			shuffle(all(t), rng);
			double tol = 1e-9 * (double)lim * (double)lim * (double)lim;
			assert(abs(signedPolyVolume(toDouble(p), t) - area2 * (double)h / 6) <= tol);
		}
	}
	cout<<"Tests passed!"<<endl;
}
