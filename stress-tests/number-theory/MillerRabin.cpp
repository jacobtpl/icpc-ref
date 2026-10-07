#include "../utilities/template.h"

#include "../../content/number-theory/MillerRabin.h"
namespace sieve {
#include "../../content/number-theory/FastEratosthenes.h"
}

ull A[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
int afactors[] = {2, 3, 5, 13, 19, 73, 193, 407521, 299210837};

const ull MR_LIM = 1ULL << 62;

// Accurate for arbitrary 64-bit numbers
ull int128_mod_mul(ull a, ull b, ull m) { return (ull)((__uint128_t)a * b % m); }
ull int128_mod_pow(ull b, ull e, ull mod) {
	ull ans = 1;
	for (; e; b = int128_mod_mul(b, b, mod), e /= 2)
		if (e & 1) ans = int128_mod_mul(ans, b, mod);
	return ans;
}
bool oldIsPrime(ull p) {
	if (p == 2) return true;
	if (p == 1 || p % 2 == 0) return false;
	ull s = p - 1;
	while (s % 2 == 0) s /= 2;
	rep(i,0,15) {
		ull a = rand() % (p - 1) + 1, tmp = s;
		ull mod = int128_mod_pow(a, tmp, p);
		while (tmp != p - 1 && mod != 1 && mod != p - 1) {
			mod = int128_mod_mul(mod, mod, p);
			tmp *= 2;
		}
		if (mod != p - 1 && tmp % 2 == 0) return false;
	}
	return true;
}

void rec(ull div, ll num, int ind, int factors) {
	if (ind == sizeof(afactors)/sizeof(*afactors)) {
		if (factors == 1) assert(isPrime(div));
		if (factors > 1) assert(!isPrime(div));
		return;
	}
	for (;;) {
		rec(div, num, ind+1, factors);
		div *= afactors[ind];
		if (num % div != 0) break;
		factors++;
	}
}

// Independent deterministic oracle: __int128 arithmetic and the first 12
// primes as bases, which is exact for all n < 3.3e24.
bool oracle(ull n) {
	if (n < 2) return false;
	for (ull p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
		if (n % p == 0) return n == p;
	}
	ull d = n - 1; int s = 0;
	while (d % 2 == 0) d /= 2, s++;
	for (ull a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
		ull x = int128_mod_pow(a, d, n);
		if (x == 1 || x == n - 1) continue;
		bool comp = true;
		rep(i,1,s) {
			x = int128_mod_mul(x, x, n);
			if (x == n - 1) { comp = false; break; }
		}
		if (comp) return false;
	}
	return true;
}

const ull DOC_LIM = 7000000000000000000ULL; // documented limit, 7e18

void check(ull n) {
	if (n > DOC_LIM) return;
	bool a = isPrime(n), b = oracle(n);
	if (a != b) {
		cout << "isPrime(" << n << ") = " << a << ", expected " << b << endl;
		exit(1);
	}
}

void testHard() {
	mt19937_64 rng(12345);
	auto rnd = [&](ull lo, ull hi) { return uniform_int_distribution<ull>(lo, hi)(rng); };
	// well-known strong pseudoprimes / Carmichael numbers
	for (ull n : {2047ULL, 1373653ULL, 9080191ULL, 25326001ULL, 3215031751ULL, 4759123141ULL,
			1122004669633ULL, 2152302898747ULL, 3474749660383ULL, 341550071728321ULL,
			3825123056546413051ULL, 561ULL, 1105ULL, 1729ULL, 2465ULL, 4681ULL, 8911ULL}) {
		assert(!isPrime(n));
		check(n);
	}
	// largest primes below 2^k
	assert(isPrime((1ULL << 61) - 1));
	assert(isPrime((1ULL << 62) - 57));
	assert(isPrime(1000000007) && isPrime(998244353) && isPrime(1000000000000000003ULL));
	// tiny numbers and both ends of the documented range
	rep(i,0,300000) check(i), check(DOC_LIM - i), check((1ULL << 62) + i - 150000);
	// Carmichael numbers (6k+1)(12k+1)(18k+1) and their non-Carmichael siblings
	for (ull k = 1; (6*k+1) * (12*k+1) <= DOC_LIM / (18*k+1); k++)
		check((6*k+1) * (12*k+1) * (18*k+1));
	// typical strong pseudoprime shapes (k+1)(r*k+1), squares, cubes, semiprimes
	rep(it,0,400000) {
		for (ull r : {2, 3, 4, 5, 7, 13}) {
			ull hi = (ull)sqrtl((long double)DOC_LIM / (long double)r) - 2;
			ull k = it % 4 == 0 ? rnd(1, 100000) : rnd(1, hi);
			check((k + 1) * (r * k + 1));
		}
		ull p = rnd(2, 2645751311ULL);
		check(p * p);
		ull q = rnd(2, 1912931);
		check(q * q * q);
		ull a = rnd(3037000000ULL / 2, 3037000000ULL) | 1, b = rnd(a, 2 * a) | 1;
		if (b <= DOC_LIM / a) check(a * b);
		// uniformly random, random bit length, and powers of two +- small
		check(rnd(0, DOC_LIM));
		check(rnd(0, DOC_LIM) >> rnd(0, 62));
		check((1ULL << rnd(1, 62)) + rnd(0, 200) - 100);
	}
	// make sure large primes are really hit: next primes after random points
	int primes = 0;
	rep(it,0,3000) {
		ull n = rnd(DOC_LIM / 2, DOC_LIM - 10000) | 1;
		while (!oracle(n)) check(n), n += 2;
		assert(isPrime(n)); primes++;
	}
	assert(primes == 3000);
}

const int MAXPR = 1e6;
int main() {
	testHard();
	auto prs = sieve::eratosthenes();
	vector<bool> isprime(MAXPR);
	for (auto i: prs) isprime[i] = true;
	for(auto &a: A) rec(1, a, 0, 0);

	rep(n,0,MAXPR) {
		if (isPrime(n) != isprime[n]) {
			cout << "fails for " << n << endl;
			return 1;
		}
	}

	ull n = 1;
	rep(i,0,1000000) {
		n ^= (ull)rand();
		n *= 1237618231ULL;
		if (n < MR_LIM && oldIsPrime(n) != isPrime(n)) {
			cout << "differs from old for " << n << endl;
			cout << "old says " << oldIsPrime(n) << endl;
			cout << "new says " << isPrime(n) << endl;
			assert(false);
		}
	}
	cout << "Tests passed!" << endl;
}
