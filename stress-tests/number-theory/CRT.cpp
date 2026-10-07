#include "../utilities/template.h"

#include "../../content/number-theory/CRT.h"

ll rmod(ll a, ll b) { return (a % b + b) % b; }

std::mt19937_64 rng(2);

ll randExp() {
	int e = rand() % 63; // 64
	return uniform_int_distribution<ll>(0, (ll)((1ULL << e) - 1))(rng);
}

int main() {
	// exhaustive on tiny inputs against brute force
	rep(m,1,25) rep(n,1,25) rep(a,-m+1,m) rep(b,-n+1,n) {
		ll l = m / __gcd(m, n) * n, exp = -1;
		rep(x,0,l) if (rmod(x - a, m) == 0 && rmod(x - b, n) == 0) {
			assert(exp == -1); // unique
			exp = x;
		}
		assert((exp != -1) == ((a - b) % __gcd(m, n) == 0));
		if (exp != -1) assert(crt(a, m, b, n) == exp);
	}
	// a, b far outside (-m, m): still a valid solution
	rep(m,1,13) rep(n,1,13) rep(a,-40,40) rep(b,-40,40)
		if ((a - b) % __gcd(m, n) == 0) {
			ll r = crt(a, m, b, n);
			assert(rmod(r - a, m) == 0 && rmod(r - b, n) == 0);
		}
	// large coprime moduli with m*n just below 2^62, known answer
	rep(it,0,1000000) {
		ll m = (1LL << 31) - (ll)(rng() % 1000), n = (1LL << 31) - (ll)(rng() % 1000);
		if (it % 3 == 0) m = (1LL << 61) - (ll)(rng() % 1000), n = 1 + (ll)(rng() % 2);
		if (it % 3 == 1) m = (1LL << 42) - (ll)(rng() % 1000), n = (1LL << 20) - (ll)(rng() % 1000);
		if (rng() % 2) swap(m, n);
		ll l = m / __gcd(m, n) * n;
		ll x = (ll)(rng() % (unsigned long long)l);
		ll a = x % m, b = x % n;
		if (rng() % 2 && a) a -= m;
		if (rng() % 2 && b) b -= n;
		assert(crt(a, m, b, n) == x);
	}
	rep(it,0,10000000) {
		ll a = randExp() * (rand() % 2 ? 1 : -1);
		ll b = randExp() * (rand() % 2 ? 1 : -1);
		ll m = randExp() + 1;
		ll n = randExp() + 1;
		ll g = __gcd(m, n);
		if (n * (__int128_t)m > LLONG_MAX) continue;
		if (n * (__int128_t)m > (1LL << 62) && (abs(a) > m || abs(b) > n)) continue;
		if ((a - b) % g == 0) {
			ll r = crt(a, m, b, n);
			if (rmod(r, m) != rmod(a, m) || rmod(r, n) != rmod(b, n)) {
				cout << a << endl;
				cout << b << endl;
				cout << m << endl;
				cout << n << endl;
			}
			assert(rmod(r, m) == rmod(a, m));
			assert(rmod(r, n) == rmod(b, n));
			if (-m < a && a < m && -n < b && b < n) {
				assert(0 <= r);
				assert(r < m*n/g);
			}
		}
	}
	cout<<"Tests passed!"<<endl;
}
