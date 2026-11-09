#include "../utilities/template.h"

#include "../../content/number-theory/ModInt.h"
#include "../../content/numerical/SubsetConvolution.h"

typedef unsigned long long ull;

template<class T> vector<T> naiveZeta(const vector<T>& a) {
	vector<T> r(sz(a));
	rep(S,0,sz(a)) rep(U,0,sz(a)) if ((U & S) == U) r[S] += a[U];
	return r;
}
template<class T> vector<T> naiveConv(const vector<T>& a, const vector<T>& b) {
	vector<T> c(sz(a));
	rep(S,0,sz(a)) for (int U = S;; U = (U - 1) & S) {
		c[S] += a[U] * b[S ^ U];
		if (!U) break;
	}
	return c;
}
template<class T> vector<T> randVec(int N, int mode) {
	vector<T> a(N);
	for (T& x : a) {
		if (mode == 0) x = T(ull(rand()) << 32 ^ ull(rand()));
		else if (mode == 1) x = T(rand() % 3);
		else if (mode == 2) x = rand() % 4 ? T(0) : T(rand());
		else x = T(ull(-1) - rand() % 3);
	}
	return a;
}
template<class T> void test(int n, int mode) {
	int N = 1 << n;
	vector<T> a = randVec<T>(N, mode), b = randVec<T>(N, mode);
	vector<T> z = a; zeta(z);
	assert(z == naiveZeta(a));
	zeta(z, 1);
	assert(z == a);
	z = a; zeta(z, 1); zeta(z);
	assert(z == a);
	vector<T> c = subsetConv(a, b);
	assert(c == naiveConv(a, b));
	assert(c == subsetConv(b, a));
}

int main() {
	rep(it,0,300) rep(n,0,9) {
		test<mint>(n, it % 4);
		test<ull>(n, it % 4);
	}

	// N = 1
	assert(subsetConv(vector<mint>{3}, vector<mint>{5}) == vector<mint>{15});
	assert(subsetConv(vector<ull>{3}, vector<ull>{5}) == vector<ull>{15});

	// all ones: c[S] = 2^|S|
	{
		int n = 12, N = 1 << n;
		vector<ull> a(N, 1), c = subsetConv(a, a);
		rep(S,0,N) assert(c[S] == 1ULL << __builtin_popcount(S));
	}

	// identity (indicator of the empty set)
	{
		int n = 10, N = 1 << n;
		vector<mint> a = randVec<mint>(N, 0), e(N); e[0] = 1;
		assert(subsetConv(a, e) == a && subsetConv(e, a) == a);
	}

	// signed answers via ull: exact when the true answer fits in ll
	{
		int n = 7, N = 1 << n;
		vector<ll> a(N), b(N);
		for (ll& x : a) x = rand() % 2000001 - 1000000;
		for (ll& x : b) x = rand() % 2000001 - 1000000;
		vector<ull> ua(all(a)), ub(all(b)), c = subsetConv(ua, ub);
		vector<ll> ref = naiveConv(a, b);
		rep(S,0,N) assert((ll)c[S] == ref[S]);
	}

	// large n: zeta and inverse round-trip
	{
		int n = 16, N = 1 << n;
		vector<ull> a = randVec<ull>(N, 0), z = a;
		zeta(z); zeta(z, 1);
		assert(z == a);
		vector<ull> ones(N, 1); zeta(ones);
		rep(S,0,N) assert(ones[S] == 1ULL << __builtin_popcount(S));
	}
	cout<<"Tests passed!"<<endl;
}
