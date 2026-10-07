#include "../utilities/template.h"

#include "../../content/data-structures/Matrix.h"

typedef unsigned long long ull;
const ll mod = 1000000007;
struct Mod {
	ll x;
	Mod(ll xx = 0) : x(xx) {}
	Mod operator+(Mod b) const { return Mod((x + b.x) % mod); }
	Mod operator*(Mod b) const { return Mod(x * b.x % mod); }
	Mod& operator+=(Mod b) { return *this = *this + b; }
	bool operator==(Mod b) const { return x == b.x; }
};

// int that aborts on overflow, to check that no product beyond the requested
// power is ever formed (signed overflow is UB for plain int / ll).
struct Chk {
	int x;
	Chk(int xx = 0) : x(xx) {}
	Chk operator*(Chk b) const {
		int r;
		assert(!__builtin_mul_overflow(x, b.x, &r));
		return Chk(r);
	}
	Chk& operator+=(Chk b) {
		assert(!__builtin_add_overflow(x, b.x, &x));
		return *this;
	}
};

mt19937_64 rng(4242);

// T = ull: wrap-around arithmetic is a ring, so results must match exactly.
template<int N> void test(int iters) {
	typedef Matrix<ull, N> M;
	auto rnd = [&]() {
		M a;
		int mode = (int)(rng() % 3);
		rep(i,0,N) rep(j,0,N)
			a.d[i][j] = mode == 0 ? rng() : mode == 1 ? rng() % 3 : -(rng() % 5);
		return a;
	};
	auto mul = [&](const M& a, const M& b) {
		M c;
		rep(i,0,N) rep(j,0,N) {
			ull s = 0;
			rep(k,0,N) s += a.d[i][k] * b.d[k][j];
			c.d[i][j] = s;
		}
		return c;
	};
	M id;
	rep(i,0,N) id.d[i][i] = 1;
	rep(it,0,iters) {
		M a = rnd(), b = rnd();
		assert((a * b).d == mul(a, b).d);
		assert((a * id).d == a.d && (id * a).d == a.d);
		vector<ull> v(N), w(N);
		for (auto& x : v) x = rng();
		rep(i,0,N) rep(j,0,N) w[i] += a.d[i][j] * v[j];
		assert(a * v == w);
		M p = id;
		rep(e,0,20) { // a^e against repeated multiplication
			assert((a ^ e).d == p.d);
			p = mul(p, a);
		}
		ll x = (ll)(rng() >> 2), y = (ll)(rng() >> 2); // huge exponents
		assert((a ^ (x + y)).d == mul(a ^ x, a ^ y).d);
		assert((a ^ LLONG_MAX).d == mul(a ^ (LLONG_MAX - 1), a).d);
	}
}

int main(int argc, char**) {
	if (argc > 1) { // benchmark mode: ./a.out bench
		auto now = [] { return chrono::steady_clock::now(); };
		auto secs = [&](auto t0) { return chrono::duration<double>(now() - t0).count(); };
		{
			static Matrix<Mod, 100> a, r;
			rep(i,0,100) rep(j,0,100) a.d[i][j] = (ll)(rng() % mod);
			auto t0 = now();
			r = a ^ (ll)1e18;
			cerr << "Matrix<Mod,100> ^ 1e18: " << secs(t0) << " s (" << r.d[0][0].x << ")\n";
		}
		{
			static Matrix<ull, 200> a, r;
			rep(i,0,200) rep(j,0,200) a.d[i][j] = rng();
			auto t0 = now();
			r = a ^ 1000;
			cerr << "Matrix<ull,200> ^ 1000: " << secs(t0) << " s (" << r.d[0][0] << ")\n";
		}
		{
			static Matrix<double, 300> a, r;
			rep(i,0,300) rep(j,0,300) a.d[i][j] = 1e-3 * (double)(rng() % 1000);
			auto t0 = now();
			rep(i,0,10) r = a * a;
			cerr << "Matrix<double,300> a*a x10: " << secs(t0) << " s (" << r.d[0][0] << ")\n";
		}
		return 0;
	}
	test<1>(2000); test<2>(5000); test<3>(5000); test<4>(2000); test<7>(500); test<16>(50);
	{ // the usage example from the header, and Fibonacci mod p
		int N = 5;
		Matrix<int, 3> A;
		A.d = {{{{1,2,3}}, {{4,5,6}}, {{7,8,9}}}};
		vector<int> vec = {1,2,3};
		vec = (A^N) * vec;
		assert((vec == vi{953856, 2160108, 3366360}));
		Matrix<Chk, 3> C; // same computation, A^5 fits in int with lots of room
		rep(i,0,3) rep(j,0,3) C.d[i][j] = A.d[i][j];
		assert(((C^N) * vector<Chk>{1,2,3})[2].x == 3366360);
		Matrix<Chk, 1> two;
		two.d[0][0] = 2;
		rep(e,0,31) assert((two ^ e).d[0][0].x == 1 << e);
		Matrix<Mod, 2> F;
		F.d = {{{{1,1}}, {{1,0}}}};
		assert((F ^ 0).d[0][0] == 1 && (F ^ 0).d[0][1] == 0);
		assert((F ^ 90).d[0][1].x == 2880067194370816120LL % mod);
		assert((F ^ (ll)1e18).d[0][1].x == 209783453);
		Matrix<double, 2> D;
		D.d = {{{{0.5, 0.5}}, {{0.25, 0.75}}}};
		auto P = D ^ 200; // stationary distribution (1/3, 2/3)
		assert(abs(P.d[0][0] - 1.0 / 3) < 1e-12 && abs(P.d[1][1] - 2.0 / 3) < 1e-12);
	}
	cout << "Tests passed!" << endl;
}
