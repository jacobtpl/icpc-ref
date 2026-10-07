#include "../utilities/template.h"

#include "../../content/number-theory/ModSum.h"

ll rmod(ll x, ll m) {
	x %= m;
	return x < 0 ? x + m : x;
}

ll rdiv(ll x, ll y) {
	x -= rmod(x, y);
	return x / y;
}

ll modsum_naive(ll to, ll c, ll k, ll m) {
	ll res = 0;
	for (int i = 0; i < (int)to; ++i)
		res += rmod(c + k * i, m);
	return res;
}

ll divsum_naive(ll to, ll c, ll k, ll m) {
	ll res = 0;
	for (int i = 0; i < (int)to; ++i)
		res += rdiv(c + k * i, m);
	return res;
}

void compare() {
	rep(to,0,30) {
		rep(c,-30,30) {
			rep(k,-30,30) {
				rep(m,1,30) {
					ll a = modsum(to, c, k, m);
					ll b = modsum_naive(to, c, k, m);
					if (a != b) {
						cout << "differ! " << to << ' ' << c << ' ' << k << ' ' << m << ": " << a << " vs " << b << endl;
						assert(false);
					}
				}
			}
		}
	}
}

void compare2() {
	rep(to,0,30) {
		rep(c,0,30) {
			rep(k,0,30) {
				rep(m,1,30) {
					ll a = divsum(to, c, k, m);
					ll b = divsum_naive(to, c, k, m);
					if (a != b) {
						cout << "differ! " << to << ' ' << c << ' ' << k << ' ' << m << ": " << a << " vs " << b << endl;
						assert(false);
					}
				}
			}
		}
	}
}

typedef unsigned __int128 u128;
typedef __int128 i128;
mt19937_64 rng(12345);
ull rnd(ull lo, ull hi) { return lo + rng() % (hi - lo + 1); }
// random value with a random bit length in [0, bits]
ull rndBits(int bits) { int b = (int)rnd(0, bits); return b ? rnd(0, (b == 64 ? 0 : 1ULL << b) - 1) : 0; }

// independent oracle (iterative floor sum): sum_{i<n} floor((a*i+b)/m) mod 2^128
u128 floorSum(u128 n, u128 m, u128 a, u128 b) {
	u128 ans = 0;
	for (;;) {
		if (a >= m) ans += n * (n - 1) / 2 * (a / m), a %= m;
		if (b >= m) ans += n * (b / m), b %= m;
		u128 y = a * n + b;
		if (y < m) break;
		n = y / m, b = y % m; swap(m, a);
	}
	return ans;
}

// large c, k, m with small `to` against brute force in 128-bit arithmetic
void bigBrute() {
	rep(it,0,200000) {
		ull to = rnd(0, 60);
		ll m = (ll)rndBits(56) + 1;
		ll c = (ll)rndBits(62) * (rng() % 2 ? 1 : -1);
		ll k = (ll)rndBits(62) * (rng() % 2 ? 1 : -1);
		i128 exp = 0;
		rep(i,0,(int)to) exp += (((i128)k * i + c) % m + m) % m;
		assert((i128)modsum(to, c, k, m) == exp);
		ull uc = rndBits(63), uk = rndBits(63), um = rndBits(57) + 1;
		u128 e2 = 0;
		rep(i,0,(int)to) e2 += ((u128)uk * i + uc) / um;
		assert(divsum(to, uc, uk, um) == (ull)e2);
		assert((ull)floorSum(to, um, uk, uc) == (ull)e2);
	}
}

// large `to` (to * m < 2^62) against the independent floor sum
void bigOracle() {
	rep(it,0,300000) {
		int tb = (int)rnd(0, 61);
		ull to = rndBits(tb), um = rndBits(61 - tb) + 1;
		ull hi = 1ULL << (61 - tb);
		if (it % 4 == 0) um = hi - min(hi - 1, rnd(0, 3)); // near the limit
		ull uc = rndBits(63), uk = rndBits(63);
		assert(divsum(to, uc, uk, um) == (ull)floorSum(to, um, uk, uc));
		ll m = (ll)um;
		ll c = (ll)rndBits(62) * (rng() % 2 ? 1 : -1);
		ll k = (ll)rndBits(62) * (rng() % 2 ? 1 : -1);
		ull c2 = (ull)((c % m + m) % m), k2 = (ull)((k % m + m) % m);
		// a, b < m here, so the 128-bit floor sum (< to^2) is exact
		i128 exp = (i128)to * c2 + (to ? (i128)k2 * ((i128)to * (to - 1) / 2) : 0)
			- (i128)m * (i128)floorSum(to, um, k2, c2);
		assert(exp >= 0 && exp <= (i128)to * (m - 1));
		assert((i128)modsum(to, c, k, m) == exp);
	}
	// edge cases
	assert(modsum(0, 5, 7, 3) == 0 && divsum(0, 5, 7, 3) == 0);
	assert(modsum(1, -1, 0, 7) == 6 && modsum(10, 3, 5, 1) == 0);
	assert(modsum(1ULL << 31, -1, 0, 1LL << 31) == ((1LL << 31) - 1) << 31);
	assert(divsum(1ULL << 32, 0, 1, 1) == ((1ULL << 32) - 1) << 31);
}

int main() {
	compare(); compare2();
	bigBrute(); bigOracle();
	assert(modsum((ll)1e18, 1, 2, 3) == (ll)1e18);
	rep(i,0,50) {
		ll t = (ll)rand() << 3;
		ll c = (ll)rand() << 2;
		ll k = (ll)rand() << 2;
		ll m = (ll)rand() >> 2;
		assert(abs(modsum(t, c, k, m) / ((long double)m/2 * t) - 1)<1e-5);
	}
	cout<<"Tests passed!"<<endl;
	return 0;
}
