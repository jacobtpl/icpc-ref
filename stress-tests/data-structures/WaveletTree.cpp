#include "../utilities/template.h"

#include "../../content/data-structures/WaveletTree.h"

template<class T>
void check(const vector<T>& a, const vector<T>& probes) {
	Wavelet<T> w(a);
	int n = sz(a);
	rep(l,0,n+1) rep(r,l,n+1) {
		vector<T> s(a.begin() + l, a.begin() + r);
		sort(all(s));
		rep(k,0,r-l) assert(w.kth(l, r, k) == s[k]);
		for (T v : probes) {
			int c = (int)(lower_bound(all(s), v) - s.begin());
			assert(w.count(l, r, v) == c);
		}
		for (T v : s) {
			int c = (int)(lower_bound(all(s), v) - s.begin());
			assert(w.count(l, r, v) == c);
		}
	}
}

template<class T>
void checkInt(mt19937& rng, int n, T lo, T hi) {
	vector<T> a(n);
	uniform_int_distribution<T> dist(lo, hi);
	for (T& x : a) x = dist(rng);
	vector<T> probes = {numeric_limits<T>::min(), numeric_limits<T>::max(),
		(T)(numeric_limits<T>::min() + 1), (T)(numeric_limits<T>::max() - 1), 0, lo, hi};
	if (n) {
		T mx = *max_element(all(a)), mn = *min_element(all(a));
		if (mx < numeric_limits<T>::max()) probes.push_back((T)(mx + 1));
		if (mn > numeric_limits<T>::min()) probes.push_back((T)(mn - 1));
	}
	check(a, probes);
}

int main() {
	mt19937 rng(2);
	// empty array
	check(vi(), vi{INT_MIN, -1, 0, 1, INT_MAX});

	rep(it,0,3000) {
		int n = (int)(rng() % 13);
		int mode = (int)(rng() % 6);
		if (mode == 0) checkInt<int>(rng, n, -3, 3);
		else if (mode == 1) checkInt<int>(rng, n, INT_MIN, INT_MIN + 4);
		else if (mode == 2) checkInt<int>(rng, n, INT_MAX - 4, INT_MAX);
		else if (mode == 3) checkInt<ll>(rng, n, LLONG_MIN, LLONG_MAX);
		else if (mode == 4) checkInt<ll>(rng, n, -5, 5);
		else checkInt<int>(rng, n, 7, 7);
	}

	// sizes around powers of two, with permutations / sorted / reversed / duplicates
	for (int n : {1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100}) {
		vi a(n);
		rep(i,0,n) a[i] = i;
		shuffle(all(a), rng);
		check(a, vi{-1, 0, n / 2, n - 1, n});
		sort(all(a));
		check(a, vi{-1, n});
		reverse(all(a));
		check(a, vi{0, n - 1});
		rep(i,0,n) a[i] = (int)(rng() % 3) * 1000000000 - 1000000000;
		check(a, vi{INT_MIN, -1000000001, -1000000000, 0, 1000000000, INT_MAX});
	}

	// non-integer types
	{
		vector<string> a = {"b", "", "abc", "a", "b", "zz", "", "ab"};
		check(a, vector<string>{"", "a", "aa", "c", "zzz"});
		vector<double> d = {0.5, -1.5, 1e18, -1e18, 0.5, 3.25, 0.0};
		check(d, vector<double>{-2e18, 0.0, 0.5, 0.75, 2e18});
	}

	// larger random test against a brute force on random queries
	{
		int n = 2000;
		vector<ll> a(n);
		for (ll& x : a) x = (ll)(rng() % 50) - 25;
		Wavelet<ll> w(a);
		rep(it,0,20000) {
			int l = (int)(rng() % (n + 1)), r = (int)(rng() % (n + 1));
			if (l > r) swap(l, r);
			vector<ll> s(a.begin() + l, a.begin() + r);
			sort(all(s));
			if (l < r) {
				int k = (int)(rng() % (r - l));
				assert(w.kth(l, r, k) == s[k]);
			}
			ll v = (ll)(rng() % 60) - 30;
			assert(w.count(l, r, v) == (int)(lower_bound(all(s), v) - s.begin()));
		}
	}
	cout<<"Tests passed!"<<endl;
}
