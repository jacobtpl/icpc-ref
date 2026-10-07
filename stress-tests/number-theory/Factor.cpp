#include "../utilities/template.h"

#include "../../content/number-theory/Factor.h"

mt19937_64 uni(12345);
void assertValid(ull N, vector<ull> prFac){
	ull cur=1;
	for (auto i: prFac){
		if (!isPrime(i)){
			cout<<N<<endl;
			cout<<i<<endl;
			assert(isPrime(i));
		}
		cur *= i;
	}
	if (cur!= N)
		cout<<cur<<' '<<N<<endl;
	assert(cur == N);
}

// independent oracle: trial division
vector<ull> naive(ull n) {
	vector<ull> r;
	for (ull p = 2; p * p <= n; p++) while (n % p == 0) r.push_back(p), n /= p;
	if (n > 1) r.push_back(n);
	return r;
}
void check(ull n) {
	auto res = factor(n);
	assertValid(n, res);
}
ull randPrime(ull lo, ull hi) {
	for (;;) {
		ull p = lo + uni() % (hi - lo + 1);
		if (isPrime(p)) return p;
	}
}

int main() {
	assert(factor(1).empty());
	assert(factor(2) == vector<ull>{2});
	assert((factor(2299) == vector<ull>{11, 19, 11}));
	// Regression: with a fixed f(x) = x^2 + 1 the cycles mod 352523 and mod
	// 352817 have the same length (821) for every starting point, so pollard()
	// never terminates on their product.
	{
		auto res = factor(352523ULL * 352817);
		sort(all(res));
		assert((res == vector<ull>{352523, 352817}));
	}
	rep(n,2,1e5) {
		auto res = factor(n);
		assertValid(n, res);
		res = factor(n*ll(n));
		assertValid(n*ll(n), res);
	}
	// compare against trial division, as multisets
	rep(n,1,300000) {
		auto res = factor(n), exp = naive(n);
		sort(all(res));
		assert(res == exp);
	}
	rep(i,0,3000) {
		ull n = 1 + uni() % (1ULL << (10 + i % 40));
		auto res = factor(n), exp = naive(n);
		sort(all(res));
		assert(res == exp);
	}
	rep(i,2,1e5) {
		ull n = 1 + (uni()%(3ul<<61));
		auto res = factor(n);
		assertValid(n, res);
	}
	rep(i,0,1e5) {
		// max number that modmul can handle
		ull n = 7268172458553106874 - i;
		auto res = factor(n);
		assertValid(n, res);
	}
	// prime powers
	rep(p,2,2000) if (isPrime(p)) {
		ull n = p;
		while (n <= 7000000000000000000ULL / p) n *= p, check(n);
	}
	// squares and cubes of larger primes
	rep(i,0,300) {
		ull p = randPrime(1e5, 2600000000ULL);
		auto res = factor(p * p);
		assert((res == vector<ull>{p, p}));
		ull q = randPrime(1e3, 1900000);
		res = factor(q * q * q);
		assert((res == vector<ull>{q, q, q}));
	}
	// semiprimes of every size, balanced and unbalanced
	rep(i,0,20000) {
		int b1 = 2 + i % 30, b2 = 2 + (i / 30) % 30;
		ull p = randPrime(1ULL << (b1 - 1), 1ULL << b1);
		ull q = randPrime(1ULL << (b2 - 1), 1ULL << b2);
		auto res = factor(p * q);
		sort(all(res));
		assert((res == vector<ull>{min(p, q), max(p, q)}));
	}
	// worst case: two primes just below sqrt(7e18)
	rep(i,0,300) {
		ull p = randPrime(2500000000ULL, 2640000000ULL);
		ull q = randPrime(2500000000ULL, 2640000000ULL);
		auto res = factor(p * q);
		sort(all(res));
		assert((res == vector<ull>{min(p, q), max(p, q)}));
	}
	// Carmichael numbers and other strong-pseudoprime-like composites
	for (ull n : {561ULL, 1105ULL, 1729ULL, 41041ULL, 825265ULL, 321197185ULL,
			5394826801ULL, 232250619601ULL, 9746347772161ULL,
			1436697831295441ULL, 60977817398996785ULL,
			3215031751ULL, 3825123056546413051ULL, 4ULL, 8ULL, 9ULL, 16ULL, 27ULL})
		check(n);
	cout<<"Tests passed!"<<endl;
}
