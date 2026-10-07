#include "../utilities/template.h"

#include "../../content/number-theory/ModMulLL.h"

const int ITERS = 5'000'000; // (not really enough to say much, need >1e10 for any kind of certainty)

ull double_modmul(ull a, ull b, ull M) {
	ll ret = a * b - M * ull(1. / (double)M * (double)a * (double)b);
	return ret + M * (ret < 0) - M * (ret >= (ll)M);
}

ull int128_modmul(ull a, ull b, ull m) { return (ull)((__uint128_t)a * b % m); }

void test(ull lim, bool expectSuccess, bool useDoubles) {
	mt19937_64 rng(1);
	uniform_int_distribution<ull> uni(1, lim);
	uniform_int_distribution<ull> uniSmall(0, lim / 10000);

	for (int i = 0;; i++) {
		if (expectSuccess && i >= ITERS) break;
		// if (i % 1'000'000 == 0) cerr << '.' << flush;
		ull c = i&1 ? lim - uniSmall(rng) : uni(rng);
		ull a = i&2 ? c - uniSmall(rng) : i&4 && !useDoubles ? (1ULL << 62) - uniSmall(rng) : uni(rng);
		ull b = i&8 ? c - uniSmall(rng) : uni(rng);
		if (a > c || b > c) continue;
		ull l = int128_modmul(a, b, c);
		ull r = useDoubles ? double_modmul(a, b, c) : modmul(a, b, c);
		if (l != r) {
			if (!expectSuccess) break;
			cout << a << ' ' << b << ' ' << c << endl;
			cout << l << ' ' << r << endl;
			abort();
		}
	}
}

ull int128_modpow(ull b, ull e, ull m) {
	ull ans = 1 % m;
	for (; e; b = int128_modmul(b, b, m), e /= 2)
		if (e & 1) ans = int128_modmul(ans, b, m);
	return ans;
}

void testEdges() {
	const ull lim = 7268172458553106874ULL, doc = 7200000000000000000ULL;
	// exhaustive for small moduli, including a = c and b = c
	rep(c,1,150) rep(a,0,c+1) rep(b,0,c+1) {
		assert(modmul(a, b, c) == ull(a * b % c));
		if (c > 1 && b < 40) {
			ull r = 1;
			rep(i,0,b) r = r * a % c;
			assert(modpow(a, b, c) == r);
		}
	}
	// extreme values around the documented bound
	mt19937_64 rng(2);
	vector<ull> cs = {1, 2, 3, (1ULL << 31) - 1, 1ULL << 32, (1ULL << 32) + 1, 1000000007,
		(1ULL << 52) - 1, 1ULL << 52, (1ULL << 53) + 1, (1ULL << 61) - 1, 1ULL << 62,
		(1ULL << 62) + 1, doc - 1, doc, lim - 1, lim};
	for (ull c : cs) {
		vector<ull> xs = {0, 1, 2, c / 2, c / 2 + 1, c - 1, c, c / 3, (ull)sqrtl((long double)c)};
		rep(i,0,2000) xs.push_back(rng() % (c + 1));
		for (ull a : xs) if (a <= c) for (ull b : {xs[0], xs[1], xs[2], xs[3], xs[4], xs[5], xs[6], xs[7], xs[8], xs[9 + a % 2000]}) {
			if (b > c) continue;
			if (modmul(a, b, c) != int128_modmul(a, b, c)) {
				cout << "modmul(" << a << "," << b << "," << c << ") = " << modmul(a, b, c) << endl;
				abort();
			}
		}
	}
	// modpow against __int128, exponents up to 2^64-1
	rep(it,0,200000) {
		ull c = it % 3 == 0 ? cs[rng() % sz(cs)] : it % 3 == 1 ? rng() % doc + 1 : (rng() >> (rng() % 62)) % doc + 1;
		ull b = it % 7 == 0 ? c : it % 7 == 1 ? c - 1 : rng() % (c + 1);
		ull e = it % 5 == 0 ? rng() % 4 : it % 5 == 1 ? ~0ULL - rng() % 4 : rng() >> (rng() % 64);
		if (c == 1 && e == 0) continue; // modpow(b, 0, 1) returns 1, not 0
		ull r = modpow(b, e, c), w = int128_modpow(b, e, c);
		if (r != w) {
			cout << "modpow(" << b << "," << e << "," << c << ") = " << r << ", expected " << w << endl;
			abort();
		}
	}
}

int main() {
	testEdges();
	const ull limDoubles = 1ULL << 52;
	test(limDoubles, true, true);
	test((ull)(limDoubles * 1.02L), false, true);

	const ull lim = 7268172458553106874ULL; // floor((sqrt(177) - 7) / 16 * 2**64)
	test(lim, true, false);
	test((ull)(lim * 1.01L), false, false);
	// test((ull)(lim * 1.001L), false);
	cout << "Tests passed!" << endl;
}
