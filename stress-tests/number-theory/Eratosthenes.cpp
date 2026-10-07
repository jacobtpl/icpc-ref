#include "../utilities/template.h"

namespace dynamic {
vi eratosthenes(int LIM) {
	const int S = (int)round(sqrt(LIM)), R = LIM / 2;
	vi pr({2}), sieve(S + 1); pr.reserve(LIM / (int)log(LIM));
	vector<array<int, 2>> cp;
	for (int i = 3; i <= S; i += 2) if (!sieve[i]) {
		cp.push_back({i, i * i / 2});
		for (int j = i * i; j <= S; j += 2 * i) sieve[j] = 1;
	}
	for (int L = 1; L <= R; L += S) {
		vector<bool> block(S);
		// array<bool, S> block{};
		for (auto &[p, idx] : cp)
			for (int i=idx; i < S+L; idx = (i+=p)) block[i-L] = 1;
		rep(i,0,min(S, R - L))
			if (!block[i]) pr.push_back((L + i) * 2 + 1);
	}
	return pr;
}
}
#include "../../content/number-theory/FastEratosthenes.h"
#include "../../content/number-theory/Eratosthenes.h"


bool naivePrime(int n) {
	if (n < 2) return 0;
	for (int d = 2; d * d <= n; d++) if (n % d == 0) return 0;
	return 1;
}

int main() {
	vi pr1 = eratosthenesSieve(LIM);
	vi pr2 = eratosthenes();
	assert(pr1 == pr2);

	for (int lim=121; lim<1000; lim++) {
		vi pr = eratosthenesSieve(lim);
		vi r = dynamic::eratosthenes(lim);
		assert(pr == r);
	}

	// Eratosthenes.h against trial division, every small lim (including
	// 0, 1, 2), checking both the returned list and isprime[0..lim).
	rep(lim,0,3000) {
		vi pr = eratosthenesSieve(lim), exp;
		rep(i,0,lim) {
			assert(isprime[i] == naivePrime(i));
			if (naivePrime(i)) exp.push_back(i);
		}
		assert(pr == exp);
	}
	// the segmented sieve (as a copy with runtime LIM) for every small LIM >= 3
	rep(lim,3,3000) {
		vi pr = dynamic::eratosthenes(lim), exp;
		rep(i,0,lim) if (naivePrime(i)) exp.push_back(i);
		assert(pr == exp);
	}
	// full documented size
	vi prAll = eratosthenesSieve(MAX_PR);
	{
		assert(sz(prAll) == 348513 && prAll.back() == 4999999);
		vector<char> comp(MAX_PR);
		for (ll i = 2; i < MAX_PR; i++) if (!comp[i])
			for (ll j = i * i; j < MAX_PR; j += i) comp[j] = 1;
		rep(i,0,MAX_PR) assert(isprime[i] == (i >= 2 && !comp[i]));
	}
	// limits that are prime, prime+1, squares of primes, powers of two, ...
	for (int lim : {4999999, 4999998, 2221*2221, 2221*2221+1, 2236*2236,
			1<<22, (1<<22)+1, 999983, 999984, 1000000}) {
		vi pr = eratosthenesSieve(lim);
		assert(pr == vi(prAll.begin(), lower_bound(all(prAll), lim)));
		int cnt = 0;
		rep(i,0,lim) cnt += isprime[i];
		assert(cnt == sz(pr));
		for (int p : pr) assert(isprime[p]);
	}
	cout<<"Tests passed!"<<endl;
}
