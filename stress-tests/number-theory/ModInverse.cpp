#include "../utilities/template.h"

ll modpow(ll a, ll e, ll mod) {
	if (e == 0) return 1;
	ll x = modpow(a * a % mod, e >> 1, mod);
	return e & 1 ? x * a % mod : x;
}
bool isPrime(ll x) {
	if (x <= 1) return false;
	for (ll i = 2; i*i <= x; ++i) {
		if (x % i == 0) return false;
	}
	return true;
}
void test(const ll mod, const int LIM) {
	assert(LIM <= mod && isPrime(mod)); // documented preconditions
	#include "../../content/number-theory/ModInverse.h"
	rep(i,1,LIM) {
		assert(0 < inv[i] && inv[i] < mod);
		assert(inv[i] * i % mod == 1);
		if (mod < 2000) assert(inv[i] == modpow(i, mod-2, mod));
	}
	delete[] (inv + 1);
}
int main() {
	// every prime below 1000 with every admissible LIM (including LIM = mod)
	rep(mod,2,1000) if (isPrime(mod))
		rep(LIM,2,mod+1) test(mod, LIM);
	// large primes, including the ones commonly used and ones close to 2^31
	mt19937 rng(7);
	for (ll mod : {1000003LL, 998244353LL, 1000000007LL, 1000000009LL, 2147483647LL, 2147483629LL, 4294967291LL}) {
		test(mod, 2);
		test(mod, 200000);
		test(mod, (int)min(mod, 2000000LL));
		rep(it,0,20) test(mod, (int)(rng() % 5000 + 2));
	}
	cout<<"Tests passed!"<<endl;
}
