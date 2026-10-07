#include "../utilities/template.h"

#include "../../content/combinatorial/NimProduct.h"

// independent oracle 1: the mex definition
const int N = 64;
int tab[N][N];
void buildMex() {
	rep(a,0,N) rep(b,0,N) {
		bitset<4 * N> s;
		rep(x,0,a) rep(y,0,b) s[tab[x][b] ^ tab[a][y] ^ tab[x][y]] = 1;
		int m = 0; while (s[m]) m++;
		tab[a][b] = m;
	}
}
// independent oracle 2: recursive Karatsuba over GF(2^(2^k))
ul rec(ul a, ul b, int bits = 64) {
	if (bits == 1) return a & b;
	int h = bits / 2;
	ul m = (1ULL << h) - 1;
	ul a0 = a & m, a1 = a >> h, b0 = b & m, b1 = b >> h;
	ul lo = rec(a0, b0, h), hi = rec(a1, b1, h);
	ul mid = rec(a0 ^ a1, b0 ^ b1, h) ^ lo;
	return (mid << h) ^ lo ^ rec(hi, 1ULL << (h - 1), h);
}

int main() {
	mt19937_64 rng(777);
	buildMex();
	rep(a,0,N) rep(b,0,N) {
		assert(ul(nb(a) * nb(b)) == (ul)tab[a][b]);
		assert(rec(a, b) == (ul)tab[a][b]);
	}
	rep(a,0,256) rep(b,0,256) {
		assert(ul(nb(a) * nb(b)) == rec(a, b));
		assert(P.x[a][b] == rec(a, b));
	}
	// powers of two, incl. 2^(2^k) * 2^(2^k) = 3/2 * 2^(2^k)
	rep(i,0,64) rep(j,0,64)
		assert(ul(nb(1ULL << i) * nb(1ULL << j)) == rec(1ULL << i, 1ULL << j));
	rep(k,0,6) {
		ul f = 1ULL << (1 << k);
		assert(ul(nb(f) * nb(f)) == (f ^ (f >> 1)));
	}
	assert(ul(nb(1ULL << 32) * nb(1ULL << 32)) == (3ULL << 31));
	auto gen = [&]() -> ul {
		ul x = rng();
		switch (rng() % 6) {
			case 0: return x >> (rng() % 64);
			case 1: return ~0ULL >> (rng() % 64);
			case 2: return 1ULL << (rng() % 64);
			case 3: return x & (x >> 7) & (x << 5);
			case 4: return ~0ULL << (rng() % 64);
			default: return x;
		}
	};
	Precalc Q = P; // mult<64> is non-const
	rep(it,0,200000) {
		ul a = gen(), b = gen(), c = gen();
		nb A(a), B(b), C(c);
		ul ab = ul(A * B);
		assert(ab == rec(a, b));
		assert(ab == ul(B * A));
		if (it % 16 == 0) assert(ab == Q.mult<64>(a, b));
		assert(ul((A * B) * C) == ul(A * (B * C)));
		assert(ul(A * (B + C)) == (ab ^ ul(A * C)));
		assert(ul(A + B) == (a ^ b));
		assert(ul(A * nb(1)) == a && ul(A * nb(0)) == 0);
		if (a) {
			nb I = inv(A);
			assert(ul(I * A) == 1);
			assert(ul(I * (A * B)) == b);
			assert(ul(pow(A, ~0ULL)) == 1); // a^(2^64-1) = 1
		}
		if (it % 8 == 0) {
			ul e = rng() % 40; nb r(1);
			rep(i,0,(int)e) r = r * A;
			assert(ul(pow(A, e)) == ul(r));
			ul e1 = gen(), e2 = gen();
			assert(ul(pow(A, e1) * pow(A, e2)) == ul(pow(pow(A, e1), 1) * pow(A, e2)));
			if (a && e1 + e2 >= e1)
				assert(ul(pow(A, e1) * pow(A, e2)) == ul(pow(A, e1 + e2)));
		}
	}
	assert(ul(inv(nb(0))) == 0 && ul(pow(nb(0), 0)) == 1);
	assert(ul(nb()) == 0);
	cout << "Tests passed!" << endl;
}
