#include "../utilities/template.h"

#include "../../content/numerical/LinearRecurrence.h"

template<class F>
void gen(vector<ll>& v, int at, F f) {
	if (at == sz(v)) f();
	else {
		rep(i,0,mod) {
			v[at] = i;
			gen(v, at+1, f);
		}
	}
}

// Oracle for huge k: power of the companion matrix.
typedef vector<vector<ll>> Mat;
Mat mul(const Mat& a, const Mat& b) {
	int n = sz(a);
	Mat c(n, vector<ll>(n));
	rep(i,0,n) rep(j,0,n) rep(k,0,n)
		c[i][j] = (c[i][j] + a[i][k] * b[k][j]) % mod;
	return c;
}
ll matRec(const vector<ll>& S, const vector<ll>& tr, ll k) {
	int n = sz(tr);
	if (n == 0) return 0;
	Mat M(n, vector<ll>(n)), R(n, vector<ll>(n));
	rep(i,0,n) R[i][i] = 1;
	rep(i,0,n-1) M[i][i+1] = 1; // (S[i..i+n-1]) -> (S[i+1..i+n])
	rep(j,0,n) M[n-1][n-1-j] = tr[j] % mod;
	for (; k; k /= 2, M = mul(M, M)) if (k & 1) R = mul(R, M);
	ll res = 0;
	rep(j,0,n) res = (res + R[0][j] * S[j]) % mod;
	return res;
}

int main() {
	mt19937_64 rng(987654321);
	auto rnd = [&](ll lo, ll hi) {
		return uniform_int_distribution<ll>(lo, hi)(rng); };
	// The order-0 recurrence S[i] = 0 (tr is empty). This is what
	// berlekampMassey returns for an all-zero sequence.
	rep(it,0,100) {
		vector<ll> S(rnd(0, 3), 0);
		ll k = it < 10 ? it : rnd(0, (ll)4e18);
		assert(linearRec(S, {}, k) == 0);
	}
	rep(n,1,5) {
		vector<ll> start(n);
		vector<ll> coef(n);
		int size = 10*n + 3;
		vector<ll> full(size);
		gen(start,0,[&]() {
			gen(coef,0,[&]() {
				for(auto &x:full) x = 0;
				rep(i,0,n) full[i] = start[i];
				rep(i,n,size) rep(j,0,n) full[i] = (full[i] + coef[j] * full[i-1 - j]) % mod;
				rep(i,0,size) {
					auto v = linearRec(start, coef, i);
					assert(v == full[i]);
				}
			});
		});
	}
	// random recurrences with larger n, S longer than n, small k
	rep(it,0,1500) {
		int n = (int)rnd(1, it % 10 ? 8 : 40), size = 3 * n + 20;
		vector<ll> S(n + rnd(0, 5)), tr(n), full(size);
		for (auto& x : tr) x = rnd(0, mod - 1);
		if (it % 3 == 0) for (auto& x : tr) x = rnd(0, 1) * rnd(0, mod - 1);
		rep(i,0,n) full[i] = rnd(0, mod - 1);
		rep(i,n,size) rep(j,0,n)
			full[i] = (full[i] + tr[j] * full[i-1-j]) % mod;
		rep(i,0,sz(S)) S[i] = full[i];
		rep(i,0,size) assert(linearRec(S, tr, i) == full[i]);
	}
	// huge k against matrix exponentiation
	rep(it,0,20000) {
		int n = (int)rnd(1, 8);
		vector<ll> S(n), tr(n);
		for (auto& x : S) x = rnd(0, mod - 1);
		for (auto& x : tr) x = rnd(0, mod - 1);
		ll k = it % 3 == 0 ? rnd(0, 1000) : it % 3 == 1 ? rnd(0, (ll)1e18)
			: LLONG_MAX - 1 - rnd(0, 1000);
		assert(linearRec(S, tr, k) == matRec(S, tr, k));
	}
	// header usage example: Fibonacci numbers (mod 5)
	{
		ll a = 0, b = 1;
		rep(k,0,200) {
			assert(linearRec({0, 1}, {1, 1}, k) == a);
			ll c = (a + b) % mod; a = b; b = c;
		}
	}
	cout<<"Tests passed!"<<endl;
}
