#include "../utilities/template.h"

#include "../../content/geometry/Point3D.h"

typedef __int128 L;
const double PI = acos(-1.0);

int main() {
	mt19937_64 rng(9);
	auto ri = [&](ll a, ll b) { return uniform_int_distribution<ll>(a, b)(rng); };
	auto rd = [&](double a, double b) { return uniform_real_distribution<double>(a, b)(rng); };
	{
		Point3D<ll> z;
		assert(z.x == 0 && z.y == 0 && z.z == 0);
		Point3D<ll> o(1, 2);
		assert(o.x == 1 && o.y == 2 && o.z == 0);
	}
	// Point3D<ll>: exact, small values (many ties) and |coord| <= 1e9
	for (ll M : {1LL, 10LL, 1000000000LL}) rep(it,0,1000000) {
		typedef Point3D<ll> P;
		ll a[3], b[3], d = ri(-M, M);
		rep(i,0,3) a[i] = ri(-M, M), b[i] = ri(-M, M);
		P p(a[0], a[1], a[2]), q(b[0], b[1], b[2]);
		assert((p < q) == (make_tuple(a[0], a[1], a[2]) < make_tuple(b[0], b[1], b[2])));
		assert((p == q) == (a[0] == b[0] && a[1] == b[1] && a[2] == b[2]));
		P s = p + q, t = p - q, m = p * d, c = p.cross(q);
		ll sv[3] = {s.x, s.y, s.z}, tv[3] = {t.x, t.y, t.z}, mv[3] = {m.x, m.y, m.z};
		ll cv[3] = {c.x, c.y, c.z};
		L dot = 0, d2 = 0;
		rep(i,0,3) {
			assert(sv[i] == a[i] + b[i] && tv[i] == a[i] - b[i]);
			assert((L)mv[i] == (L)a[i] * d);
			int j = (i + 1) % 3, k = (i + 2) % 3;
			assert((L)cv[i] == (L)a[j] * b[k] - (L)a[k] * b[j]);
			dot += (L)a[i] * b[i], d2 += (L)a[i] * a[i];
		}
		if (d) {
			P e = p / d;
			assert(e.x == a[0] / d && e.y == a[1] / d && e.z == a[2] / d);
		}
		assert((L)p.dot(q) == dot && (L)p.dist2() == d2);
		assert(abs(p.dist() - sqrt((double)d2)) <= 1e-9 * (double)M);
		// aliasing (operators take const references)
		P x = p; x = x + x; assert(x == p * 2);
		x = p; x = x - x; assert(x == P());
		x = p; x = x.cross(x); assert(x == P());
		assert(p.dot(p) == p.dist2());
	}
	// Point3D<double>
	for (double M : {1.0, 1e3, 1e6}) rep(it,0,300000) {
		typedef Point3D<double> P;
		P v(rd(-M, M), rd(-M, M), rd(-M, M)), ax(rd(-M, M), rd(-M, M), rd(-M, M));
		double len = v.dist(), eps = 1e-9 * M;
		assert(abs(len - sqrt(v.x * v.x + v.y * v.y + v.z * v.z)) <= eps);
		if (len < 1e-3 * M || ax.dist() < 1e-3 * M) continue;
		P u = v.unit();
		assert(abs(u.dist() - 1) < 1e-12 && (u * len - v).dist() <= eps);
		// spherical coordinates
		double ph = v.phi(), th = v.theta();
		assert(-PI <= ph && ph <= PI && 0 <= th && th <= PI);
		P w(len * sin(th) * cos(ph), len * sin(th) * sin(ph), len * cos(th));
		assert((w - v).dist() <= eps);
		// normal
		P n = v.normal(ax);
		if (v.cross(ax).dist() > 1e-3 * M * M) {
			assert(abs(n.dist() - 1) < 1e-9 && abs(n.dot(v)) <= eps && abs(n.dot(ax)) <= eps);
			assert(n.dot(v.cross(ax)) > 0);
		}
		// rotation against an explicit rotation matrix
		double a = rd(-10, 10), s = sin(a), c = cos(a);
		P k = ax / ax.dist();
		double kk[3] = {k.x, k.y, k.z}, in[3] = {v.x, v.y, v.z}, out[3];
		double K[3][3] = {{0, -kk[2], kk[1]}, {kk[2], 0, -kk[0]}, {-kk[1], kk[0], 0}};
		rep(i,0,3) {
			out[i] = 0;
			rep(j,0,3) {
				double k2 = 0; // (K^2)[i][j]
				rep(l,0,3) k2 += K[i][l] * K[l][j];
				out[i] += ((i == j) + s * K[i][j] + (1 - c) * k2) * in[j];
			}
		}
		P r = v.rotate(a, ax);
		assert((r - P(out[0], out[1], out[2])).dist() <= eps);
		assert(abs(r.dist() - len) <= eps && abs(r.dot(k) - v.dot(k)) <= eps);
		assert((r.rotate(-a, ax) - v).dist() <= eps);
		assert((ax.rotate(a, ax) - ax).dist() <= eps);
		// around the z axis this is a ccw rotation of (x, y)
		P rz = v.rotate(a, P(0, 0, 5));
		assert(abs(rz.x - (v.x * c - v.y * s)) <= eps && abs(rz.y - (v.x * s + v.y * c)) <= eps
			&& abs(rz.z - v.z) <= eps);
	}
	{
		typedef Point3D<double> P;
		assert(P(0, 0, 1).theta() == 0 && abs(P(0, 0, -1).theta() - PI) < 1e-15);
		assert(abs(P(1, 0, 0).theta() - PI / 2) < 1e-15 && P(1, 0, 0).phi() == 0);
		assert(abs(P(0, -1, 0).phi() + PI / 2) < 1e-15 && P(-1, 0, 0).phi() == PI);
		P r = P(1, 0, 0).rotate(PI / 2, P(0, 0, 1));
		assert((r - P(0, 1, 0)).dist() < 1e-15);
		assert((P(1, 0, 0).normal(P(0, 1, 0)) - P(0, 0, 1)).dist() < 1e-15);
	}
	cout<<"Tests passed!"<<endl;
}
