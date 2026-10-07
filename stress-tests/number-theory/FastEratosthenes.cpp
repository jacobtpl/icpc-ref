#include "../utilities/template.h"

#include "../../content/number-theory/FastEratosthenes.h"

bool naivePrime(int n) {
	if (n < 2) return 0;
	for (int d = 2; d * d <= n; d++) if (n % d == 0) return 0;
	return 1;
}

int main() {
	assert(LIM == 1000000);
	vi pr = eratosthenes();
	// independent oracle: plain byte sieve
	vector<char> comp(LIM);
	vi exp;
	for (ll i = 2; i < LIM; i++) if (!comp[i]) {
		exp.push_back((int)i);
		for (ll j = i * i; j < LIM; j += i) comp[j] = 1;
	}
	assert(sz(exp) == 78498 && exp.back() == 999983);
	assert(pr == exp);
	rep(i,0,LIM) assert(isPrime[i] == (i >= 2 && !comp[i]));
	rep(i,0,20000) assert(isPrime[i] == naivePrime(i));
	rep(i,LIM-20000,LIM) assert(isPrime[i] == naivePrime(i));
	// calling it again must give the same result
	assert(eratosthenes() == exp);
	rep(i,0,LIM) assert(isPrime[i] == (i >= 2 && !comp[i]));
	cout<<"Tests passed!"<<endl;
}
