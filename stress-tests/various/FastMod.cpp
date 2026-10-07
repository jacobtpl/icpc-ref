#include "../utilities/template.h"

#include "../../content/various/FastMod.h"

typedef unsigned long long ull;
struct OldBarrett {
	ull b, m;
	OldBarrett(ull b) : b(b), m(-1ULL / b) {}
	ull reduce(ull a) {
		ull q = (ull)((__uint128_t(m) * a) >> 64), r = a - q * b;
		return r >= b ? r - b : r;
	}
};

// If EIGHT is defined, we compute eight simultaneous factorials, thus measuring
// throughput instead of latency.
// #define EIGHT

#ifdef EIGHT

// Compute k * (p-1)! % mod p. Should equal -k.
#define TEST() \
	INIT(0) \
	INIT(1) \
	INIT(2) \
	INIT(3) \
	INIT(4) \
	INIT(5) \
	INIT(6) \
	INIT(7) \
	for (int i = 1; i < mod; i++) { \
		UPDATE(0) \
		UPDATE(1) \
		UPDATE(2) \
		UPDATE(3) \
		UPDATE(4) \
		UPDATE(5) \
		UPDATE(6) \
		UPDATE(7) \
	} \
	FINISH(0) \
	FINISH(1) \
	FINISH(2) \
	FINISH(3) \
	FINISH(4) \
	FINISH(5) \
	FINISH(6) \
	FINISH(7)

#else

#define TEST() \
	INIT(0) \
	for (int i = 1; i < mod; i++) { \
		UPDATE(0) \
	} \
	FINISH(0)

#endif

void perf_plain(int mod) {
#define INIT(x) ll ret##x = (x + 1);
#define UPDATE(x) ret##x = (ret##x * i) % mod;
#define FINISH(x) cout << ret##x << endl;
	TEST()
#undef INIT
#undef UPDATE
#undef FINISH
}

template<int mod>
void perf_const() {
#define INIT(x) ll ret##x = (x + 1);
#define UPDATE(x) ret##x = (ret##x * i) % mod;
#define FINISH(x) cout << ret##x << endl;
	TEST()
#undef INIT
#undef UPDATE
#undef FINISH
}

#define INIT(x) ll ret##x = (x + 1);
#define UPDATE(x) ret##x = bar.reduce(ret##x * i);
#define FINISH(x) cout << ret##x << endl;
void perf_old_barrett(int mod) {
	OldBarrett bar(mod);
	TEST()
}

void perf_barrett(int mod) {
	FastMod bar(mod);
	TEST()
}
#undef INIT
#undef UPDATE
#undef FINISH

ull rand_u64() {
	ull ret = rand();
	ret <<= 23;
	ret ^= rand();
	ret <<= 23;
	ret ^= rand();
	return ret;
}

#ifndef BENCH
#define main1 main
#endif

// Correctness
int main1() {
	const int bflim = 3000;
	rep(a,0,bflim) rep(b,2,bflim) {
		FastMod bar(b);
		ull ret = bar.reduce(a);
		assert((ret == 0) == (a == 0));
		if (ret >= (ull)b) ret -= b;
		assert(ret == (ull)(a % b));
	}
	// Edge cases: b = 1, powers of two, b and a near 2^64, a around multiples of b.
	// The result must be congruent to a mod b, lie in [0, 2b), and be 0 only if a = 0.
	{
		vector<ull> bs = {1, 2, 3, 5, 7, 10, 1000000007, 998244353, (1ULL << 32) - 1, 1ULL << 32, (1ULL << 32) + 1,
			(1ULL << 63) - 1, 1ULL << 63, (1ULL << 63) + 1, -1ULL, -2ULL, -3ULL, 6700417, 4294967291ULL};
		rep(i,1,64) bs.push_back(1ULL << i), bs.push_back((1ULL << i) - 1), bs.push_back((1ULL << i) + 1);
		rep(it,0,2000) bs.push_back(rand_u64() >> (rand() % 64)), bs.push_back(-1ULL - (rand_u64() >> (rand() % 64)));
		for (ull b : bs) {
			if (b == 0) continue;
			FastMod bar(b);
			vector<ull> as = {0, 1, 2, b - 1, b, b + 1, 2 * b - 1, 2 * b, 2 * b + 1, -1ULL, -2ULL, 1ULL << 63,
				(1ULL << 63) - 1, -1ULL / b * b, -1ULL / b * b - 1, -1ULL / b * b + 1, b * b, b * b - 1, (b - 1) * (b - 1)};
			rep(it,0,50) as.push_back(rand_u64()), as.push_back(rand_u64() / b * b), as.push_back(rand_u64() / b * b - 1),
				as.push_back(rand_u64() >> (rand() % 64));
			for (ull a : as) {
				ull ret = bar.reduce(a);
				assert((ret == 0) == (a == 0));
				assert((__uint128_t)ret < (__uint128_t)2 * b);
				assert(ret % b == a % b);
			}
		}
		// typical usage: modular multiplication of reduced values
		for (ull b : {1000000007ULL, 998244353ULL, 2ULL, 1ULL, 4294967295ULL}) {
			FastMod bar(b);
			rep(it,0,200000) {
				ull x = rand_u64() % b, y = rand_u64() % b, r = bar.reduce(x * y);
				assert(r < 2 * b && r % b == x * y % b);
			}
		}
	}
	rep(it,0,10'000'000) {
		ull a = rand_u64();
		ull b = rand_u64();
		if (b == 0) continue;
		FastMod bar(b);
		ull ret = bar.reduce(a);
		if (ret >= b) ret -= b;
		assert(ret == a % b);
	}
	cout<<"Tests passed!"<<endl;
	return 0;
}

#ifndef MOD
#define MOD 90217093 //202171241
#endif
int mod = MOD;

// Performance; compile with -DBENCH. Runs all four variants, or the given one.
int main2(int argc, char** argv) {
	rep(which,0,4) {
		if (argc > 1 && atoi(argv[1]) != which) continue;
		auto t0 = chrono::steady_clock::now();
		if (which == 0) perf_plain(mod); // 7.529 for 8, 1.714 for 1
		if (which == 1) perf_const<MOD>(); // 0.971 for 8, 0.499 for 1
		if (which == 2) perf_old_barrett(mod); // 1.094 for 8, 0.564 for 1
		if (which == 3) perf_barrett(mod); // 0.870 for 8, 0.405 for 1
		auto t1 = chrono::steady_clock::now();
		const char* names[] = {"plain %", "compile-time constant %", "old barrett", "FastMod"};
		cerr << names[which] << ": " << chrono::duration<double>(t1 - t0).count() << " s" << endl;
	}
	return 0;
}

#ifdef BENCH
int main(int argc, char** argv) { return main2(argc, argv); }
#endif
