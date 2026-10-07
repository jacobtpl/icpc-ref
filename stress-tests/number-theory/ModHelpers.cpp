#include "../utilities/template.h"

#include "../../content/number-theory/ModHelpers.h"

ll mpow(ll b, ll e) {
	ll r = 1;
	for (b %= MOD; e; e /= 2, b = b * b % MOD)
		if (e & 1) r = r * b % MOD;
	return r;
}

int main() {
	assert(MV == 2000010 && MOD == 998244353);
	init();
	// inv / fact / ifact over the whole table
	ll f = 1;
	assert(fact[0].v == 1 && ifact[0].v == 1);
	rep(i,1,MV) {
		f = f * i % MOD;
		assert(0 < inv[i].v && inv[i].v < MOD);
		assert((ll)inv[i].v * i % MOD == 1);
		assert(fact[i].v == f);
		assert((ll)ifact[i].v * f % MOD == 1);
	}
	rep(i,1,3000) assert(inv[i].v == mpow(i, MOD - 2));
	// Pascal's triangle
	const int N = 600;
	vector<vector<ll>> C(N, vector<ll>(N));
	rep(n,0,N) rep(k,0,n+1) {
		C[n][k] = k == 0 || k == n ? 1 : (C[n-1][k-1] + C[n-1][k]) % MOD;
		assert(choose(n, k).v == C[n][k]);
	}
	// large arguments against an independently computed product formula
	mt19937 rng(5);
	rep(it,0,300000) {
		int n = it % 3 == 0 ? MV - 1 - (int)(rng() % 50) : (int)(rng() % MV);
		int k = it % 5 == 0 ? (int)(rng() % 4) : it % 5 == 1 ? n - (int)(rng() % min(n + 1, 4)) : (int)(rng() % (n + 1));
		ll want;
		if (min(k, n - k) <= 30) {
			ll num = 1, den = 1;
			rep(i,0,min(k, n - k)) num = num * (n - i) % MOD, den = den * (i + 1) % MOD;
			want = num * mpow(den, MOD - 2) % MOD;
		} else {
			// C(n,k) = C(n-1,k-1) * n / k
			want = (ll)choose(n - 1, k - 1).v * n % MOD * mpow(k, MOD - 2) % MOD;
		}
		assert(choose(n, k).v == want);
		assert(choose(n, k).v == choose(n, n - k).v);
	}
	assert(choose(0, 0).v == 1 && choose(MV - 1, 0).v == 1 && choose(MV - 1, MV - 1).v == 1);
	assert(choose(MV - 1, 1).v == MV - 1);
	// init() is idempotent
	init();
	assert(choose(10, 3).v == 120 && fact[20].v == 2432902008176640000LL % MOD);
	cout<<"Tests passed!"<<endl;
}
