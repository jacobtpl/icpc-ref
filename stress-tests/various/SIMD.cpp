#include "../utilities/template.h"

#include "../../content/various/SIMD.h"

mt19937 rng(4242);
int rnd(int lo, int hi) { return uniform_int_distribution<int>(lo, hi)(rng); }

ll slowDot(int n, short* a, short* b) {
	ll r = 0;
	rep(i,0,n) if (a[i] < b[i]) r += a[i] * b[i];
	return r;
}

int main() {
	if (!__builtin_cpu_supports("avx2")) {
		cout << "Tests passed! (skipped: no AVX2)" << endl;
		return 0;
	}
	// helpers
	assert(all_zero(zero()) && !all_one(zero()));
	assert(all_one(one()) && !all_zero(one()));
	assert(sumi32(zero()) == 0 && sumi32(one()) == -8);
	rep(it,0,200000) {
		int v[8]; ll s = 0; bool az = 1, ao = 1;
		int mode = rnd(0, 5);
		rep(i,0,8) {
			v[i] = mode == 0 ? 0 : mode == 1 ? -1 : mode == 2 ? (rnd(0, 7) ? 0 : 1 << rnd(0, 31))
				: mode == 3 ? ~(rnd(0, 7) ? 0 : 1 << rnd(0, 31)) : rnd(-1000000, 1000000);
			s += v[i]; az &= v[i] == 0; ao &= v[i] == -1;
		}
		mi m = L(v[0]);
		if (INT_MIN <= s && s <= INT_MAX) assert(sumi32(m) == s);
		assert(all_zero(m) == az);
		assert(all_one(m) == ao);
	}
	// example_filteredDotProduct against the scalar loop
	rep(it,0,300000) {
		int n = it < 200000 ? rnd(0, 70) : rnd(0, 600);
		int mode = rnd(0, 4);
		vector<short> a(n + 1), b(n + 1);
		auto gen = [&]() -> short {
			if (mode == 0) return (short)rnd(-32768, 32767);
			if (mode == 1) return (short)rnd(0, 32767);
			if (mode == 2) return (short)rnd(-32768, -1);
			if (mode == 3) return (short)(rnd(0, 1) ? rnd(32760, 32767) : rnd(-32768, -32760));
			return (short)rnd(-3, 3);
		};
		rep(i,0,n) a[i] = gen(), b[i] = gen();
		int o = rnd(0, 1); // unaligned start
		int m = max(0, n - o);
		ll x = example_filteredDotProduct(m, a.data() + o, b.data() + o);
		ll y = slowDot(m, a.data() + o, b.data() + o);
		if (x != y) {
			cout << "mismatch n=" << m << " mode=" << mode << " got " << x << " expected " << y << endl;
			return 1;
		}
	}
	// extremes
	{
		int n = 1 << 20;
		vector<short> a(n, -32768), b(n, -32767);
		assert(example_filteredDotProduct(n, a.data(), b.data()) == slowDot(n, a.data(), b.data()));
		a.assign(n, 32766), b.assign(n, 32767);
		assert(example_filteredDotProduct(n, a.data(), b.data()) == slowDot(n, a.data(), b.data()));
		a.assign(n, -32768), b.assign(n, 32767);
		assert(example_filteredDotProduct(n, a.data(), b.data()) == slowDot(n, a.data(), b.data()));
	}
	cout << "Tests passed!" << endl;
}
