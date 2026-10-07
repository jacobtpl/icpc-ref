#include "../utilities/template.h"

#include "../../content/numerical/MatrixInverse-mod.h"

typedef vector<vector<ll>> vvll;
mt19937_64 rng(20240601);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }
ll norm(ll x) { return (x % mod + mod) % mod; }

// Oracle: rank by plain row elimination with partial pivoting.
int rankOracle(vvll a) {
	int n = sz(a), r = 0;
	for (auto& v : a) for (auto& x : v) x = norm(x);
	rep(c,0,n) {
		int p = r;
		while (p < n && !a[p][c]) p++;
		if (p == n) continue;
		swap(a[r], a[p]);
		ll inv = modpow(a[r][c], mod - 2);
		rep(j,r+1,n) {
			ll f = a[j][c] * inv % mod;
			rep(k,c,n) a[j][k] = norm(a[j][k] - f * a[r][k]);
		}
		r++;
	}
	return r;
}

void fail(const vvll& m, const char* what) {
	cerr << what << endl;
	for (auto& r : m) { for (ll x : r) cerr << x << ' '; cerr << endl; }
	abort();
}

void check(const vvll& m) {
	int n = sz(m);
	vvll inv = m;
	int rk = matInv(inv);
	if (rk != rankOracle(m)) fail(m, "wrong rank");
	if (rk < n) return;
	rep(i,0,n) rep(j,0,n) {
		if (inv[i][j] < 0 || inv[i][j] >= mod) fail(m, "not in [0, mod)");
		ll s1 = 0, s2 = 0; // both A * inv and inv * A must be I
		rep(k,0,n) {
			s1 = (s1 + norm(m[i][k]) * inv[k][j]) % mod;
			s2 = (s2 + inv[i][k] * norm(m[k][j])) % mod;
		}
		if (s1 != (i == j) || s2 != (i == j)) fail(m, "not an inverse");
	}
}

template<class F>
void rec(int i, int j, vvll& A, int lo, int hi, F f) {
	if (i == sz(A)) f();
	else if (j == sz(A)) rec(i+1, 0, A, lo, hi, f);
	else rep(v,lo,hi+1) {
		A[i][j] = v;
		rec(i, j+1, A, lo, hi, f);
	}
}

vvll randMat(int n, ll lo, ll hi) {
	vvll m(n, vector<ll>(n));
	for (auto& r : m) for (auto& x : r) x = rnd(lo, hi);
	return m;
}

int main() {
	// exhaustive tiny matrices: n <= 2 over {0..4}, n = 3 over {0, 1, 2},
	// n = 4 over {0, 1}
	rep(n,0,5) {
		vvll m(n, vector<ll>(n));
		rec(0, 0, m, 0, n <= 2 ? 4 : n == 3 ? 2 : 1, [&]() { check(m); });
	}
	// 1 x 1 edge values
	for (ll v : {0LL, 1LL, 2LL, (ll)mod - 1, (ll)mod - 2, (ll)(mod + 1) / 2}) {
		vvll m{{v}}, inv = m;
		int rk = matInv(inv);
		assert(rk == (v != 0));
		if (rk) assert(inv[0][0] * v % mod == 1);
	}
	rep(it,0,60000) {
		int n = (int)rnd(1, it % 100 ? 6 : 30);
		vvll m;
		switch (it % 6) {
		case 0: m = randMat(n, 0, mod - 1); break; // generic, full range
		case 1: m = randMat(n, 0, 1); break; // often singular
		case 2: // near the modulus, or negative residues in (-mod, mod)
			m = it % 12 == 2 ? randMat(n, mod - 3, mod - 1)
				: randMat(n, 1 - mod, mod - 1);
			break;
		case 3: { // sparse: many zero pivots, forces row and column swaps
			m = randMat(n, 0, mod - 1);
			for (auto& r : m) for (auto& x : r) if (rnd(0, 2)) x = 0;
			break;
		}
		case 4: { // rank r product
			int r = (int)rnd(0, n);
			vvll A = randMat(max(n, r), 0, mod - 1), B = randMat(max(n, r), 0, mod - 1);
			m.assign(n, vector<ll>(n));
			rep(i,0,n) rep(j,0,n) rep(k,0,r)
				m[i][j] = (m[i][j] + A[i][k] * B[k][j]) % mod;
			break;
		}
		default: { // permutation matrix with random nonzero scalings
			vi p(n); iota(all(p), 0); shuffle(all(p), rng);
			m.assign(n, vector<ll>(n));
			rep(i,0,n) m[i][p[i]] = rnd(1, mod - 1);
		}
		}
		check(m);
	}
#ifdef MATINV_UNREDUCED
	// Opt-in: entries must already be reduced to (-mod, mod). An entry
	// equal to mod is picked as a "nonzero" pivot and the invertible
	// matrix {{mod, 1}, {1, 0}} == {{0, 1}, {1, 0}} is reported singular.
	{
		vvll m{{mod, 1}, {1, 0}};
		int rk = matInv(m);
		cerr << "rank of {{mod, 1}, {1, 0}}: " << rk << ", want 2" << endl;
		assert(rk == 2);
	}
#endif
	// O(n^3) at n = 150
	check(randMat(150, 0, mod - 1));
	cout<<"Tests passed!"<<endl;
}
