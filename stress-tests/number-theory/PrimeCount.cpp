#include "../utilities/template.h"

#include "../../content/number-theory/PrimeCount.h"

const int LIM = 2000000;
vector<ll> cnt(LIM + 1), sum(LIM + 1);

void sieve() {
	vector<bool> comp(LIM + 1);
	rep(i,2,LIM+1) {
		if (!comp[i]) for (ll j = (ll)i * i; j <= LIM; j += i) comp[j] = 1;
		cnt[i] = cnt[i-1] + !comp[i];
		sum[i] = sum[i-1] + (comp[i] ? 0 : i);
	}
}

ll naive(ll x, bool k) { return k ? sum[x] : cnt[x]; }

void check(ll n, bool full) {
	rep(k,0,2) {
		PrimeCount pc(n, k);
		int s = pc.s;
		assert((ll)s * s <= n && (ll)(s + 1) * (s + 1) > n);
		assert(sz(pc.lo) == s + 1 && sz(pc.hi) == s + 1);
		rep(x,1,s+1) assert(pc.lo[x] == naive(x, k));
		rep(i,1,s+1) assert(pc.hi[i] == naive(n / i, k));
		assert(pc.get(n) == naive(n, k));
		for (ll i = 1; i <= n; i++) {
			ll x = n / i;
			assert(pc.get(x) == naive(x, k));
			if (!full && x <= s) break;
			i = n / x;
		}
	}
}

int main() {
	sieve();
	rep(n,1,3001) check(n, true);
	rep(k,1,1415) {
		ll sq = (ll)k * k;
		for (ll n : {sq - 1, sq, sq + 1, sq + k, sq + k - 1})
			if (n >= 1 && n <= LIM) check(n, false);
	}
	rep(it,0,1500) check(rand() % LIM + 1, false);

	PrimeCount a((ll)1e9);
	assert(a.hi[1] == 50847534);
	PrimeCount b((ll)1e9, 1);
	assert(b.hi[1] == 24739512092254535LL);
	PrimeCount c((1LL << 32) - 1);
	assert(c.hi[1] == 203280221 && c.get(c.n / 1000) == 302499);
	PrimeCount d((1LL << 32) - 1, 1);
	assert(d.hi[1] == 425649736193687430LL);
	PrimeCount e((ll)1e10);
	assert(e.hi[1] == 455052511 && e.get(1) == 0 && e.get(2) == 1);

	cout << "Tests passed!" << endl;
	return 0;
}
