#include "../utilities/template.h"

#include "../../content/number-theory/phiFunction.h"

int main() {
	calculatePhi();
	assert(phi[0] == 0 && phi[1] == 1 && phi[2] == 1);
	// definition, by brute force
	rep(n,1,3000) {
		int c = 0;
		rep(k,1,n+1) c += __gcd(n, k) == 1;
		assert(phi[n] == c);
	}
	// independent oracle for the whole table: smallest prime factor sieve
	// and phi(p^k * m) = (p-1) p^(k-1) phi(m)
	vi spf(LIM), ref(LIM);
	for (int i = 2; i < LIM; i++) if (!spf[i])
		for (int j = i; j < LIM; j += i) if (!spf[j]) spf[j] = i;
	ref[1] = 1;
	rep(n,2,LIM) {
		int p = spf[n], m = n / p;
		ref[n] = m % p ? ref[m] * (p - 1) : ref[m] * p;
	}
	rep(n,0,LIM) assert(phi[n] == ref[n]);
	// sum_{d|n} phi(d) = n
	vector<ll> s(200000);
	rep(d,1,sz(s)) for (int j = d; j < sz(s); j += d) s[j] += phi[d];
	rep(n,1,sz(s)) assert(s[n] == n);
	// calling it again must give the same table
	calculatePhi();
	rep(n,0,LIM) assert(phi[n] == ref[n]);
	cout<<"Tests passed!"<<endl;
}
