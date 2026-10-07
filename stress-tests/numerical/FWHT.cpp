#include "../utilities/template.h"

#include "../../content/numerical/FWHT.h"

mt19937_64 rng(31337);
ll rnd(ll lo, ll hi) { return lo + (ll)(rng() % (unsigned long long)(hi - lo + 1)); }

// Definition: H[i] = sum_j (-1)^popcount(i & j) a[j]
template<class T> vector<T> naive(const vector<T>& a) {
	int n = sz(a);
	vector<T> r(n);
	rep(i,0,n) rep(j,0,n)
		r[i] += __builtin_popcount(i & j) & 1 ? -a[j] : a[j];
	return r;
}

template<class T> void test(int iters, int maxk, ll lim) {
	rep(it,0,iters) {
		int k = (int)rnd(0, maxk), n = 1 << k;
		vector<T> a(n), b(n);
		for (T& x : a) x = (T)rnd(-lim, lim);
		for (T& x : b) x = (T)rnd(-lim, lim);
		auto fa = a, fb = b;
		fwht(fa); fwht(fb);
		assert(fa == naive(a));
		// inverse as documented: transform again and divide by the size
		auto back = fa;
		fwht(back);
		for (T& x : back) x /= (T)n;
		assert(back == a);
		// xor convolution
		vector<T> want(n), c(n);
		rep(i,0,n) rep(j,0,n) want[i ^ j] += a[i] * b[j];
		rep(i,0,n) c[i] = fa[i] * fb[i];
		fwht(c);
		for (T& x : c) x /= (T)n;
		assert(c == want);
	}
}

int main() {
	test<int>(20000, 6, 50);
	test<ll>(20000, 6, 1000000);
	test<double>(5000, 6, 1000); // integer valued, so exact
	test<int>(30, 10, 20);
	test<ll>(30, 10, 30000);
	{ // n = 1 is the identity
		vector<ll> one{-7};
		fwht(one);
		assert(one[0] == -7);
	}
	{ // large: involution up to the factor n, and Parseval, at n = 2^20
		int n = 1 << 20;
		vector<ll> a(n);
		for (ll& x : a) x = rnd(-1000000, 1000000);
		auto b = a;
		fwht(b);
		__int128 s1 = 0, s2 = 0;
		rep(i,0,n) s1 += (__int128)a[i] * a[i], s2 += (__int128)b[i] * b[i];
		assert(s2 == s1 * n);
		fwht(b);
		rep(i,0,n) assert(b[i] == a[i] * n);
	}
#ifdef FWHT_EMPTY
	// Opt-in: an empty vector calls __builtin_ctz(0) / __builtin_clz(0),
	// which is undefined (aborts or runs off the vector).
	{
		vector<int> e;
		fwht(e);
		assert(e.empty());
	}
#endif
	cout<<"Tests passed!"<<endl;
}
