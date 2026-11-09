#include "../utilities/template.h"

#include "../../content/number-theory/PrimitiveRoot.h"

typedef __int128_t lll;
mt19937_64 rng(5);
const ull LIM = 7000000000000000000ULL;

ull pw(ull b, ull e, ull m) {
	ull r = 1 % m;
	for (; e; e /= 2, b = (ull)((lll)b * b % m))
		if (e & 1) r = (ull)((lll)r * b % m);
	return r;
}

// multiplicative order of a mod m by repeated multiplication
ull orderBrute(ull a, ull m) {
	ull x = a % m, k = 1;
	while (x != 1) { x = x * a % m; k++; }
	return k;
}

ull smallestRootBrute(ull m, ull phi) {
	for (ull g = 1;; g++) if (__gcd(g, m) == 1 && orderBrute(g, m) == phi) return g;
}

// independent trial-division factorization (distinct primes)
vector<ull> trialFactor(ull n) {
	vector<ull> r;
	for (ull d = 2; d * d <= n; d++) if (n % d == 0) {
		r.push_back(d);
		while (n % d == 0) n /= d;
	}
	if (n > 1) r.push_back(n);
	return r;
}

bool isRoot(ull g, ull p, const vector<ull>& qs) {
	for (ull q : qs) if (pw(g, (p-1)/q, p) == 1) return 0;
	return 1;
}

// check g == smallest primitive root of prime p, given distinct prime factors qs of p-1
void checkKnown(ull p, const vector<ull>& qs) {
	ull g = primRoot(p);
	if (p == 2) { assert(g == 1); return; }
	if (!isRoot(g, p, qs) || g >= p) { cout << p << ' ' << g << endl; assert(0); }
	for (ull c = 2; c < g; c++) if (isRoot(c, p, qs)) { cout << p << ' ' << g << ' ' << c << endl; assert(0); }
}

ull randPrime(ull lo, ull hi) { // isPrime from MillerRabin.h is tested separately
	for (;;) {
		ull x = uniform_int_distribution<ull>(lo, hi)(rng);
		if (isPrime(x)) return x;
	}
}

int main() {
	assert(primRoot(2) == 1);
	assert(primRoot(3) == 2);
	assert(primRoot(5) == 2);
	assert(primRoot(7) == 3);
	assert(primRoot(65537) == 3);
	assert(primRoot(998244353) == 3);
	assert(primRoot(1000000007) == 5);
	assert(primRoot((1ULL << 61) - 1) == 37);
	assert(primRoot(40487) == 5);

	// all primes below N vs brute-force order computation
	const int N = 5000;
	vector<bool> comp(N);
	rep(p,2,N) if (!comp[p]) {
		for (int j = 2*p; j < N; j += p) comp[j] = 1;
		ull g = primRoot(p);
		assert(g == smallestRootBrute(p, p-1));
		// p^k and 2p^k statement from the description
		if (p > 2) {
			ull h = pw(g, p-1, (ull)p*p) == 1 ? g + p : g;
			for (ull pk = p; pk < 300000; pk *= p) {
				ull phi = pk / p * (p-1);
				assert(__gcd(h, pk) == 1 && orderBrute(h, pk) == phi);
				ull hh = h % 2 ? h : h + pk;
				assert(orderBrute(hh % (2*pk), 2*pk) == phi);
				if (pk > 300000 / (ull)p) break;
			}
		}
	}

	// all primes below 200000 vs trial-division factoring + __int128 powmod
	rep(p,2,200000) if (isPrime(p)) checkKnown(p, trialFactor(p-1));

	// random primes of every bit length up to 2^40 (trial division of p-1 is cheap)
	rep(b,2,41) rep(it,0,40) {
		ull p = randPrime(1ULL<<(b-1), (1ULL<<b)-1);
		checkKnown(p, trialFactor(p-1));
	}

	// large primes < 7e18 built so that the factorization of p-1 is known
	for (int cnt = 0; cnt < 100;) { // p = 2ab+1, a,b ~ 1e9
		ull a = randPrime(1e8, 1.8e9), b = randPrime(1e8, 1.8e9), p = 2*a*b + 1;
		if (p > LIM || !isPrime(p)) continue;
		checkKnown(p, {2, a, b}); cnt++;
	}
	for (int cnt = 0; cnt < 100;) { // p-1 = 2^k * m, many repeated factors
		int k = (int)(rng() % 50) + 10;
		ull m = (rng() % ((1ULL << (62 - k)) - 1)) * 2 + 1, p = (m << k) + 1;
		if (p > LIM || !isPrime(p)) continue;
		auto qs = trialFactor(m); qs.push_back(2);
		checkKnown(p, qs); cnt++;
	}
	for (int cnt = 0; cnt < 100;) { // p-1 = 2 * three primes ~1e6
		ull a = randPrime(1e5, 2e6), b = randPrime(1e5, 2e6), c = randPrime(1e5, 2e6), p = 2*a*b*c + 1;
		if (p > LIM || !isPrime(p)) continue;
		checkKnown(p, {2, a, b, c}); cnt++;
	}
	for (int cnt = 0; cnt < 30;) { // safe primes 2q+1
		ull q = randPrime(1e18, 3.4e18), p = 2*q + 1;
		if (!isPrime(p)) continue;
		checkKnown(p, {2, q}); cnt++;
	}
	// largest primes below 7e18 (the documented limit) and 2^62
	for (ull p : {6999999999999999973ULL, 4611686018427387847ULL}) {
		assert(isPrime(p));
		auto f = factor(p - 1); // Factor.h is stress-tested separately
		sort(all(f)); f.erase(unique(all(f)), f.end());
		for (ull q : f) assert(isPrime(q) && (p-1) % q == 0);
		checkKnown(p, f);
	}
	cout << "Tests passed!" << endl;
}
