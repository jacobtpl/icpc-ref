#include "../utilities/template.h"

#include "../../content/various/ConstantIntervals.h"

// Compile with -DBENCH for timings.
mt19937 rng(777);
int rnd() { return (int)(rng() >> 1); }

ll calls;
// Runs constantIntervals on v (indexed from `from`) and checks the result against a linear scan.
template<class T>
void check(int from, const vector<T>& v) {
	int n = sz(v), to = from + n;
	vector<tuple<int, int, T>> got, want;
	calls = 0;
	constantIntervals(from, to, [&](int x) {
		assert(from <= x && x < to); // never evaluates f outside the range
		calls++;
		return v[x - from];
	}, [&](int lo, int hi, T val) { got.emplace_back(lo, hi, val); });
	for (int i = 0, j; i < n; i = j) {
		for (j = i; j < n && v[j] == v[i];) j++;
		want.emplace_back(from + i, from + j, v[i]);
	}
	assert(got == want);
	// O(k log n) evaluations of f
	int k = sz(want), lg = 1;
	while ((1 << lg) < n + 1) lg++;
	assert(calls <= 2 + (ll)k * (lg + 1));
}

int main() {
	// Exhaustive: all non-decreasing sequences over {0,1,2,3} of length <= 12 (as bitmask of steps).
	rep(n,0,13) rep(mask,0,1 << max(n - 1, 0)) {
		vi v(n);
		rep(i,1,n) v[i] = min(3, v[i-1] + ((mask >> (i-1)) & 1));
		// cap at 3 values merges some masks, harmless
		check(0, v);
		check(-5, v);
		if (n == 0) break;
	}
	// Random monotone sequences, both directions, with many duplicates, various offsets and types.
	rep(it,0,200000) {
		int n = rnd() % 40, steps = rnd() % 3 ? rnd() % 4 : rnd() % 50;
		vi v(n);
		int cur = rnd() % 100 - 50;
		rep(i,0,n) {
			if (steps && rnd() % (n / steps + 1) == 0) cur += rnd() % 3 + 1;
			v[i] = cur;
		}
		if (rnd() % 2) reverse(all(v));
		int from = rnd() % 3 == 0 ? 0 : rnd() % 2001 - 1000;
		check(from, v);
		if (it % 4 == 0) {
			vector<ll> w(all(v));
			for (auto& x : w) x *= (ll)1e16;
			check(from, w);
			vector<string> s;
			for (int x : v) s.push_back(string(3, (char)('a' + (x + 60) % 26)) + to_string(x));
			check(from, s);
		}
	}
	// Edge cases: empty and reversed ranges call nothing.
	{
		int cnt = 0;
		constantIntervals(5, 5, [&](int) { cnt++; return 0; }, [&](int, int, int) { cnt++; });
		constantIntervals(7, 3, [&](int) { cnt++; return 0; }, [&](int, int, int) { cnt++; });
		assert(cnt == 0);
	}
	// Large implicit ranges (no array): floor(x / d) and a step function near INT_MAX.
	for (int n : {1, 2, 1000, 1 << 20, 1000000007, INT_MAX / 2}) for (int d : {1, 7, 1000, 123456789}) {
		if (n / d > 200000) continue;
		ll evals = 0; int expectLo = 0, k = 0;
		constantIntervals(0, n, [&](int x) { evals++; return x / d; }, [&](int lo, int hi, int val) {
			assert(lo == expectLo && lo == val * d && hi == min(n, lo + d));
			expectLo = hi; k++;
		});
		assert(expectLo == n && k == (n - 1) / d + 1);
		assert(evals <= 2 + (ll)k * 32);
	}

#ifdef BENCH
	{
		// k intervals in a range of size n: count evaluations of f and time.
		for (int n : {1000000, 1000000000}) for (int k : {1, 10, 1000, 100000, 1000000}) {
			if (k > n) continue;
			int d = n / k;
			ll evals = 0, cnt = 0;
			auto t0 = chrono::steady_clock::now();
			constantIntervals(0, n, [&](int x) { evals++; return x / d; }, [&](int, int, int) { cnt++; });
			auto t1 = chrono::steady_clock::now();
			cerr << "n=" << n << " k=" << cnt << ": " << evals << " evals of f ("
				<< (double)evals / (double)cnt << " per interval), "
				<< chrono::duration<double>(t1 - t0).count() << " s" << endl;
		}
		// worst case k = n over a vector
		int n = 5000000; vi v(n); iota(all(v), 0);
		ll evals = 0, cnt = 0;
		auto t0 = chrono::steady_clock::now();
		constantIntervals(0, n, [&](int x) { evals++; return v[x]; }, [&](int, int, int) { cnt++; });
		auto t1 = chrono::steady_clock::now();
		cerr << "all distinct n=k=5e6: " << evals << " evals, " << chrono::duration<double>(t1 - t0).count() << " s" << endl;
	}
#endif
	cout<<"Tests passed!"<<endl;
}
