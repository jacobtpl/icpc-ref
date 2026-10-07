#include "../utilities/template.h"

#include "../../content/numerical/MatrixInverse.h"

typedef vector<vector<double>> vvd;
mt19937_64 rng(555);
double rnd(double lo, double hi) {
	return uniform_real_distribution<double>(lo, hi)(rng);
}
int rndi(int lo, int hi) { return uniform_int_distribution<int>(lo, hi)(rng); }

// Oracle: exact rank of an integer matrix by elimination modulo a
// large prime (correct with overwhelming probability; all matrices
// below additionally have a rank known by construction).
int rankOracle(const vector<vector<ll>>& m) {
	const ll p = 998244353;
	int n = sz(m), r = 0;
	vector<vector<ll>> a = m;
	for (auto& v : a) for (auto& x : v) x = (x % p + p) % p;
	rep(c,0,n) {
		int q = r;
		while (q < n && !a[q][c]) q++;
		if (q == n) continue;
		swap(a[r], a[q]);
		ll inv = 1, b = a[r][c];
		for (ll e = p - 2; e; e /= 2, b = b * b % p) if (e & 1) inv = inv * b % p;
		rep(j,r+1,n) {
			ll f = a[j][c] * inv % p;
			rep(k,c,n) a[j][k] = ((a[j][k] - f * a[r][k]) % p + p) % p;
		}
		r++;
	}
	return r;
}

double residual(const vvd& A, const vvd& B) { // max |A*B - I|
	int n = sz(A); double worst = 0;
	rep(i,0,n) rep(j,0,n) {
		long double s = 0;
		rep(k,0,n) s += (long double)A[i][k] * B[k][j];
		worst = max(worst, (double)fabsl(s - (i == j)));
	}
	return worst;
}

void fail(const vvd& m, const char* what, double x) {
	cerr << what << ' ' << x << endl;
	cerr << setprecision(17);
	for (auto& r : m) { for (double v : r) cerr << v << ' '; cerr << endl; }
	abort();
}

void checkInv(const vvd& m, double tol) {
	vvd inv = m;
	int rk = matInv(inv);
	if (rk != sz(m)) fail(m, "reported singular, rank", rk);
	double r1 = residual(m, inv), r2 = residual(inv, m);
	if (!(max(r1, r2) <= tol)) fail(m, "residual", max(r1, r2));
}

int main() {
	// n = 0 and n = 1
	{
		vvd e; assert(matInv(e) == 0);
		vvd z{{0}}; assert(matInv(z) == 0);
		vvd o{{-4}}; assert(matInv(o) == 1 && o[0][0] == -0.25);
	}
	// exhaustive 0/1/-1 matrices for n <= 3: rank and inverse
	rep(n,1,4) {
		int cells = n * n, total = 1;
		rep(i,0,cells) total *= 3;
		rep(mask,0,total) {
			vvd m(n, vector<double>(n));
			vector<vector<ll>> mi(n, vector<ll>(n));
			int x = mask;
			rep(i,0,n) rep(j,0,n) mi[i][j] = x % 3 - 1, x /= 3,
				m[i][j] = (double)mi[i][j];
			vvd inv = m;
			int rk = matInv(inv);
			if (rk != rankOracle(mi)) fail(m, "wrong rank", rk);
			if (rk == n && residual(m, inv) > 1e-12) fail(m, "residual", 0);
		}
	}
	rep(it,0,40000) {
		int n = rndi(1, it % 50 ? 8 : 40);
		vvd m(n, vector<double>(n));
		vector<vector<ll>> mi(n, vector<ll>(n));
		int type = it % 5;
		if (type == 0) { // well conditioned: random plus a dominant diagonal
			rep(i,0,n) rep(j,0,n) m[i][j] = rnd(-1, 1) + (i == j) * 2.0 * n;
			checkInv(m, 1e-9);
		} else if (type == 1) { // permuted, scaled diagonal: needs pivoting
			vi p(n); iota(all(p), 0); shuffle(all(p), rng);
			rep(i,0,n) m[i][p[i]] = (rndi(0, 1) ? 1 : -1) * pow(10, rnd(-3, 3));
			checkInv(m, 1e-9);
		} else if (type == 2) { // unimodular integer matrix (det = +-1)
			rep(i,0,n) m[i][i] = 1;
			rep(s,0,2*n) {
				int a = rndi(0, n-1), b = rndi(0, n-1);
				if (a == b) continue;
				double c = rndi(-1, 1);
				if (rndi(0, 1)) rep(k,0,n) m[a][k] += c * m[b][k];
				else rep(k,0,n) m[k][a] += c * m[k][b];
			}
			if (n <= 8) checkInv(m, 1e-6);
		} else { // small integer entries, possibly rank deficient
			int r = type == 3 ? n : rndi(0, n);
			vector<vector<ll>> A(n, vector<ll>(r)), B(r, vector<ll>(n));
			for (auto& v : A) for (auto& x : v) x = rndi(-3, 3);
			for (auto& v : B) for (auto& x : v) x = rndi(-3, 3);
			if (type == 3) rep(i,0,n) rep(j,0,n) mi[i][j] = rndi(-5, 5);
			else rep(i,0,n) rep(j,0,n) rep(k,0,r) mi[i][j] += A[i][k] * B[k][j];
			rep(i,0,n) rep(j,0,n) m[i][j] = (double)mi[i][j];
			if (n > 8) continue; // keep the elimination error far below 1e-12
			vvd inv = m;
			int rk = matInv(inv), want = rankOracle(mi);
			if (rk != want) fail(m, "wrong rank, want", want);
			if (rk == n) {
				double res = residual(m, inv);
				if (res > 1e-6) fail(m, "residual", res);
			}
		}
	}
	// Hilbert matrices up to n = 6 (condition number ~1.5e7)
	rep(n,1,7) {
		vvd m(n, vector<double>(n));
		rep(i,0,n) rep(j,0,n) m[i][j] = 1.0 / (i + j + 1);
		checkInv(m, 1e-7);
	}
	// O(n^3) at n = 200
	{
		int n = 200;
		vvd m(n, vector<double>(n));
		rep(i,0,n) rep(j,0,n) m[i][j] = rnd(-1, 1);
		checkInv(m, 1e-7);
	}
	cout<<"Tests passed!"<<endl;
}
