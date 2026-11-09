#include "../utilities/template.h"

#include "../../content/data-structures/UnionFind.h"
#include "../../content/various/ParallelBinarySearch.h"

mt19937_64 rng(1337);
ll rnd(ll lo, ll hi) { return (ll)(rng() % (unsigned long long)(hi - lo + 1)) + lo; }

int rounds(int n) { int r = 0; for (int x = n + 1; x; x /= 2) r++; return r; }

// 1. Counter predicate with known thresholds; also checks the hook contract.
void testCounter(int n, int q) {
	vi th(q);
	rep(j,0,q) th[j] = (int)rnd(0, n + 1);
	int cnt = 0, resets = 0, checks = 0, applies = 0;
	auto R = [&]() { cnt = 0; resets++; };
	auto A = [&](int i) { assert(i == cnt); cnt++; applies++; };
	auto C = [&](int j) { checks++; return cnt >= th[j]; };
	vi r = parBinSearch(n, q, R, A, C);
	assert(sz(r) == q);
	assert(r == th);
	int R0 = rounds(n);
	assert(resets == R0);
	assert(applies == R0 * n);
	assert(checks <= R0 * q);
}

// 2. DSU connectivity, exactly the Usage snippet, vs naive prefix closure.
void testDSU(int n, int m, int q) {
	vi a(m), b(m), u(q), v(q);
	rep(i,0,m) a[i] = (int)rnd(0, n - 1), b[i] = (int)rnd(0, n - 1);
	rep(j,0,q) u[j] = (int)rnd(0, n - 1), v[j] = (int)rnd(0, n - 1);
	UF uf(0); // edges a[i]-b[i], queries u[j],v[j]
	auto R = [&]() { uf = UF(n); };
	auto A = [&](int i) { uf.join(a[i], b[i]); };
	auto C = [&](int j) { return uf.sameSet(u[j], v[j]); };
	vi r = parBinSearch(m, q, R, A, C);
	assert(sz(r) == q);
	rep(j,0,q) {
		int ans = m + 1;
		rep(k,0,m+1) {
			UF t(n);
			rep(i,0,k) t.join(a[i], b[i]);
			if (t.sameSet(u[j], v[j])) { ans = k; break; }
		}
		assert(r[j] == ans);
	}
}

// 3. Prefix sum >= threshold, with values near the ll range.
void testPrefix(int n, int q, bool big) {
	ll M = big ? (ll)3e18 / max(n, 1) : 3;
	vector<ll> w(n), th(q);
	rep(i,0,n) w[i] = rnd(0, M);
	rep(j,0,q) th[j] = rnd(-M, M * (n + 1));
	ll s = 0;
	auto R = [&]() { s = 0; };
	auto A = [&](int i) { s += w[i]; };
	auto C = [&](int j) { return s >= th[j]; };
	vi r = parBinSearch(n, q, R, A, C);
	assert(sz(r) == q);
	rep(j,0,q) {
		int ans = n + 1; ll p = 0;
		rep(k,0,n+1) {
			if (p >= th[j]) { ans = k; break; }
			if (k < n) p += w[k];
		}
		assert(r[j] == ans);
	}
}

int main() {
	rep(n,0,4) rep(q,0,4) { testCounter(n, q); testPrefix(n, q, 0); testPrefix(n, q, 1); }
	rep(n,1,4) rep(m,0,4) rep(q,0,4) testDSU(n, m, q);
	rep(it,0,20000) {
		int n = (int)rnd(0, 20), q = (int)rnd(0, 20);
		testCounter(n, q);
		testPrefix(n, q, it & 1);
		testDSU((int)rnd(1, 8), n, q);
	}
	rep(it,0,30) {
		int n = (int)rnd(1000, 3000), q = (int)rnd(0, 500);
		testCounter(n, q);
		testPrefix(n, q, it & 1);
	}
	testCounter(1 << 20, 1 << 17);
	testCounter((1 << 20) - 2, 1000);
	cout << "Tests passed!" << endl;
}
