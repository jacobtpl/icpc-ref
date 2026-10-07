#include "../utilities/template.h"

#include "../../content/numerical/IntDeterminant.h"

typedef vector<vector<ll>> vvll;

ll idet(vvll& a) { // integer determinant
	int n = sz(a); ll ans = 1;
	rep(i,0,n) {
		rep(j,i+1,n) {
			while (a[j][i] != 0) { // gcd step
				ll t = a[i][i] / a[j][i]; // can take mod-inv if mod p
				rep(k,i,n) a[i][k] -= a[j][k] * t;
				swap(a[i], a[j]);
				ans *= -1;
			}
		}
		if (!a[i][i]) return 0;
		ans *= a[i][i];
	}
	return ans;
}

double det(vector<vector<double>>& a) {
	int n = sz(a); double res = 1;
	rep(i,0,n) {
		int b = i;
		rep(j,i+1,n) if (fabs(a[j][i]) > fabs(a[b][i])) b = j;
		if (i != b) swap(a[i], a[b]), res *= -1;
		res *= a[i][i];
		if (res == 0) return 0;
		rep(j,i+1,n) {
			double v = a[j][i] / a[i][i];
			if (v != 0) rep(k,i+1,n) a[j][k] -= v * a[i][k];
		}
	}
	return res;
}

// Oracle 1: Leibniz formula over all permutations, exact in __int128.
__int128 leibniz(const vvll& a) {
	int n = sz(a);
	vi p(n); iota(all(p), 0);
	__int128 res = 0;
	do {
		__int128 t = 1;
		rep(i,0,n) t *= a[i][p[i]];
		rep(i,0,n) rep(j,0,i) if (p[j] > p[i]) t = -t;
		res += t;
	} while (next_permutation(all(p)));
	return res;
}

// Oracle 2: determinant modulo a prime p by field Gaussian elimination.
ll powmod(ll b, ll e, ll p) {
	ll r = 1;
	for (b %= p; e; e /= 2, b = b * b % p) if (e & 1) r = r * b % p;
	return r;
}
ll detPrime(vvll a, ll p) {
	int n = sz(a); ll res = 1;
	for (auto& r : a) for (auto& x : r) x = (x % p + p) % p;
	rep(i,0,n) {
		int b = i;
		while (b < n && !a[b][i]) b++;
		if (b == n) return 0;
		if (b != i) swap(a[i], a[b]), res = (p - res) % p;
		res = res * a[i][i] % p;
		ll inv = powmod(a[i][i], p - 2, p);
		rep(j,i+1,n) {
			ll f = a[j][i] * inv % p;
			rep(k,i,n) a[j][k] = ((a[j][k] - f * a[i][k]) % p + p) % p;
		}
	}
	return res;
}

template<class F>
void rec(int i, int j, vvll& A, F f) {
	if (i == sz(A)) {
		f();
	}
	else if (j == sz(A[i])) {
		rec(i+1, 0, A, f);
	}
	else {
		rep(v,-3,4) {
			A[i][j] = v;
			rec(i, j+1, A, f);
		}
	}
}

mt19937_64 rng(31337);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }

void fail(const vvll& m, ll got, ll want) {
	for (auto& r : m) { for (ll x : r) cerr << x << ' '; cerr << endl; }
	cerr << "got " << got << " want " << want << endl;
	abort();
}

// compare against the three prime factors of mod = 12345 = 3 * 5 * 823
void checkCRT(const vvll& m) {
	vvll m2 = m;
	ll got = det(m2);
	if (got < 0 || got >= mod) fail(m, got, -1);
	for (ll p : {3, 5, 823}) {
		ll want = detPrime(m, p);
		if (got % p != want) fail(m, got % p, want);
	}
}

int main() {
	assert(mod == 3 * 5 * 823);
	// exhaustive: all matrices with entries in [-3, 3] for n <= 3
	rep(n,0,4) {
		vvll mat(n, vector<ll>(n, 0)), mat2;
		vector<vector<double>> mat3(n, vector<double>(n, 0));
		rec(0,0,mat,[&]() {
			rep(i,0,n) rep(j,0,n) mat3[i][j] = (double)mat[i][j];
			mat2 = mat; ll h = det(mat2);
			ll a = (ll)round(det(mat3)) % mod;
			mat2 = mat; ll b = idet(mat2) % mod;
			if (a < 0) a += mod;
			if (b < 0) b += mod;
			if (a != b || h != a) fail(mat, h, b);
		});
	}
	// random small matrices against the exact Leibniz expansion
	rep(it,0,200000) {
		int n = (int)rnd(0, 6);
		ll lim = it % 4 == 0 ? 2 : it % 4 == 1 ? mod - 1 : 1000000000;
		vvll m(n, vector<ll>(n));
		for (auto& r : m) for (auto& x : r)
			x = it % 8 < 4 ? rnd(0, lim) : rnd(-lim, lim);
		if (n >= 2 && it % 5 == 0) m[rnd(0, n-1)] = m[rnd(0, n-1)];
		if (n >= 2 && it % 7 == 0) {
			int c = (int)rnd(0, n-1);
			rep(i,0,n) m[i][c] = 0;
		}
		vvll red = m; // reduced copy keeps the exact oracle within int128
		for (auto& r : red) for (auto& x : r) x %= mod;
		ll want = (ll)(leibniz(red) % mod);
		if (want < 0) want += mod;
		vvll m2 = m;
		ll got = det(m2);
		if (got != want) fail(m, got, want);
	}
	// larger matrices, structured and random, against the CRT oracle
	rep(it,0,3000) {
		int n = (int)rnd(1, 40);
		vvll m(n, vector<ll>(n));
		int type = it % 6;
		ll lim = type == 0 ? 1 : type == 1 ? mod - 1 : 1000000000;
		for (auto& r : m) for (auto& x : r)
			x = it % 2 ? rnd(0, lim) : rnd(-lim, lim);
		if (type == 3) { // low rank: product of n x r and r x n
			int r = (int)rnd(0, n);
			vvll A(n, vector<ll>(r)), B(r, vector<ll>(n));
			for (auto& v : A) for (auto& x : v) x = rnd(0, mod - 1);
			for (auto& v : B) for (auto& x : v) x = rnd(0, mod - 1);
			rep(i,0,n) rep(j,0,n) {
				m[i][j] = 0;
				rep(k,0,r) m[i][j] = (m[i][j] + A[i][k] * B[k][j]) % mod;
			}
		}
		if (type == 4) // entries are multiples of the factors of mod
			for (auto& r : m) for (auto& x : r)
				x = rnd(0, 20) * (ll[]){3, 5, 823, 15, 2469, 4115, 1}[rnd(0, 6)] % mod;
		if (type == 5) // triangular with a permutation of the rows
			{ rep(i,0,n) rep(j,0,i) m[i][j] = 0; shuffle(all(m), rng); }
		checkCRT(m);
	}
	// identity, zero and n = 0
	rep(n,0,50) {
		vvll id(n, vector<ll>(n)), z = id;
		rep(i,0,n) id[i][i] = 1;
		assert(det(id) == 1);
		assert(det(z) == (n == 0));
	}
	// O(N^3): a 400 x 400 matrix, verified with the CRT oracle
	{
		int n = 400;
		vvll m(n, vector<ll>(n));
		for (auto& r : m) for (auto& x : r) x = rnd(0, mod - 1);
		checkCRT(m);
	}
	cout<<"Tests passed!"<<endl;
}
