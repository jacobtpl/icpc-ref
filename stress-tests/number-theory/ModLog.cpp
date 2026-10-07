#include "../utilities/template.h"

#include "../../content/number-theory/ModLog.h"

ll mpow(ll b, ll e, ll m) {
	__int128 r = 1 % m, x = b % m;
	for (; e; e /= 2, x = x * x % m)
		if (e & 1) r = r * x % m;
	return (ll)r;
}

void fail(ll a, ll b, ll m, ll want, ll res) {
	cerr << "FAIL" << endl;
	cerr << "Expected log(" << a << ", " << b << ", " << m << ") = " << want << ", found " << res << endl;
	exit(1);
}

// smallest x > 0 with a^x = b (mod m) for every b, by walking the orbit of a
vector<ll> brute(ll a, int m) {
	vector<ll> ans(m, -1);
	ll b = a % m;
	rep(x,1,max(m,2)+1) {
		if (ans[b] != -1) break;
		ans[b] = x;
		b = b * a % m;
	}
	return ans;
}

int main() {
	// exhaustive for small moduli
	const int lim = 130;
	rep(m,1,lim) rep(a,0,m) {
		vector<ll> ans = brute(a, m);
		rep(b,0,m) {
			ll res = modLog(a, b, m);
			if (ans[b] != res) fail(a, b, m, ans[b], res);
		}
	}
	mt19937_64 rng(11);
	// medium moduli (highly composite, prime powers, primes): every b that is
	// reachable plus random unreachable ones
	vector<int> ms = {1024, 59049, 65536, 46656, 30030, 510510, 720720, 999983, 1000000,
		1048576, 531441, 823543, 999999, 1594323, 2097152, 1000003};
	rep(it,0,60) ms.push_back((int)(rng() % 300000) + 2);
	for (int m : ms) rep(it,0,6) {
		ll a = it == 0 ? 2 : it == 1 ? 6 : it == 2 ? m - 1 : (ll)(rng() % m);
		vector<ll> ans = brute(a, m);
		rep(q,0,300) {
			ll b = q < 150 ? mpow(a, (ll)(rng() % (2 * m)) + 1, m) : (ll)(rng() % m);
			if (q < 40) b = mpow(a, q + 1, m); // short tails for non-coprime a
			ll res = modLog(a, b, m);
			if (ans[b] != res) fail(a, b, m, ans[b], res);
		}
	}
	// Large moduli up to the largest m for which m*m fits in a signed 64-bit
	// integer. Plant a solution x and check the result is a valid, not larger,
	// solution. For prime m with a of known order q also check minimality and -1.
	const ll MLIM = 3037000499LL;
	rep(it,0,400) {
		ll m = it % 4 == 0 ? MLIM - (ll)(rng() % 1000) : it % 4 == 1 ? (ll)(rng() % MLIM) + 1
			: it % 4 == 2 ? (1LL << 31) - (ll)(rng() % 3) : (ll)(rng() % 1000000000) + 1;
		ll a = it % 8 < 2 ? (ll)(rng() % 50) : (ll)(rng() % m);
		if (it % 16 == 9) { // force many shared factors
			m = 1LL << (20 + rng() % 11); a = 2 * (ll)(rng() % (m / 2));
		}
		a %= m;
		ll x = it % 3 == 0 ? (ll)(rng() % 60) + 1 : (ll)(rng() % m) + 1;
		ll b = mpow(a, x, m);
		ll res = modLog(a, b, m);
		if (res < 1 || res > x || mpow(a, res, m) != b) fail(a, b, m, x, res);
	}
	// prime modulus p = 2q+1 with q prime: an element a = g^2 != 1 has order
	// exactly q, so the answer is x mod q (or q), and non-residues give -1.
	for (ll p : {3037000427LL, 2147483579LL, 1000000007LL, 2000000579LL}) {
		ll q = (p - 1) / 2;
		rep(i,2,60000) assert(p % i && q % i); // p, q < 3.1e9 are prime
		rep(it,0,40) {
			ll g = (ll)(rng() % (p - 3)) + 2, a = mpow(g, 2, p);
			if (a == 1) continue;
			ll x = (ll)(rng() % q) + 1, b = mpow(a, x, p);
			ll res = modLog(a, b, p);
			if (res != x) fail(a, b, p, x, res);
			ll nb = p - b; // -1 is a non-residue as p = 3 mod 4
			if ((res = modLog(a, nb, p)) != -1) fail(a, nb, p, -1, res);
			if ((res = modLog(a, 0, p)) != -1) fail(a, 0, p, -1, res);
		}
	}
#ifdef BIG_M
	// Opt-in: e * a % m overflows a signed 64-bit integer once m exceeds
	// about 3.04e9, although m is an ll and O(sqrt m) is cheap far beyond that.
	{
		ll m = 1000000000039LL, a = 2, x = 123456789, b = mpow(a, x, m);
		ll res = modLog(a, b, m);
		if (res < 1 || res > x || mpow(a, res, m) != b) fail(a, b, m, x, res);
	}
#endif
	cout<<"Tests passed!"<<endl;
}
