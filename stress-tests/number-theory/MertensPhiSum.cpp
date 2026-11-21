#include "../utilities/template.h"

#include "../../content/number-theory/MertensPhiSum.h"

const int N = 12000000; // > LIM so M/Phi recurse past the sieve
ll muSum[N], phiSum[N];

string str(lll x) {
	bool neg = x < 0; if (neg) x = -x;
	string s; do s += char('0' + (int)(x % 10)), x /= 10; while (x);
	if (neg) s += '-';
	reverse(all(s)); return s;
}

// Independent oracle: smallest-prime-factor sieve, mu and phi from it.
void naiveSieve() {
	vector<int> spf(N), m(N), ph(N);
	for (int i = 2; i < N; i++) if (!spf[i])
		for (int j = i; j < N; j += i) if (!spf[j]) spf[j] = i;
	m[1] = ph[1] = 1;
	rep(i,2,N) {
		int p = spf[i], j = i / p;
		m[i] = j % p == 0 ? 0 : -m[j];
		ph[i] = j % p == 0 ? ph[j] * p : ph[j] * (p - 1);
	}
	rep(i,1,N) muSum[i] = muSum[i-1] + m[i], phiSum[i] = phiSum[i-1] + ph[i];
}

// Independent Du recursion for Phi, with a different base: Phi(n) = T(n) - sum_{d>=2} Phi(n/d).
map<ll, lll> memo2;
lll phiDu(ll n) {
	if (n < N) return phiSum[n];
	auto it = memo2.find(n);
	if (it != memo2.end()) return it->second;
	lll r = (lll)n * (n + 1) / 2;
	for (ll l = 2, h; l <= n; l = h + 1)
		h = n / (n / l), r -= (lll)(h - l + 1) * phiDu(n / l);
	return memo2[n] = r;
}

void check(ll n) {
	ll m = M(n); lll p = Phi(n);
	if (n < N) {
		if (m != muSum[n] || p != phiSum[n]) {
			cout << "n=" << n << " M=" << m << " expected " << muSum[n]
				<< " Phi=" << str(p) << " expected " << phiSum[n] << endl;
			assert(false);
		}
	} else {
		lll q = phiDu(n);
		if (p != q) {
			cout << "n=" << n << " Phi=" << str(p) << " expected " << str(q) << endl;
			assert(false);
		}
	}
}

int main() {
	naiveSieve();
	initMu();
	assert(M(0) == 0 && M(1) == 1 && M(2) == 0 && M(3) == -1);
	assert(Phi(0) == 0 && Phi(1) == 1 && Phi(2) == 2 && Phi(3) == 4);

	// Call order should not matter: Phi before M, big before small.
	assert(Phi((ll)1e10) == phiDu((ll)1e10));
	assert(M((ll)1e10) == -33722);

	rep(n,0,300000) check(n);
	rep(d,-3,4) check(LIM + d), check(2 * LIM + d), check(N - 1 + min(d, 0));
	rep(it,0,5000) check(rand() % N);
	rep(it,0,5000) check(LIM + rand() % (N - LIM));

	// OEIS A084237 (Mertens at 10^k) and A002088 (sum of phi at 10^k).
	ll mert[] = {1, -1, 1, 2, -23, -48, 212, 1037, 1928, -222, -33722};
	const char* ph[] = {"1", "32", "3044", "304192", "30397486", "3039650754",
		"303963552392", "30396356427242", "3039635516365908",
		"303963551173008414", "30396355092886216366"};
	ll n = 1;
	rep(k,0,11) {
		assert(M(n) == mert[k]);
		assert(str(Phi(n)) == ph[k]);
		if (k >= 7) assert(Phi(n) == phiDu(n));
		n *= 10;
	}
	rep(it,0,5) check((ll)1e8 + (ll)rand() * 7 % (ll)9e8);
	assert(Phi((ll)6e9) > (lll)LLONG_MAX); // overflows ll, hence __int128
	cout << "Tests passed!" << endl;
}
