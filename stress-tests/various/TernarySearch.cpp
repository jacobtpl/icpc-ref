#include "../utilities/template.h"

#include "../../content/various/TernarySearch.h"

// Variant described in the header: non-strict on the left side,
// strict on the right, i.e. f(a) <= ... <= f(i) > ... > f(b);
// returns the largest maximizer. "Reverse the loop at (B)" has to be
// read as mirroring it (scan downwards from b, accumulating in b).
// Merely reversing the iteration order of (B) is wrong, e.g. on
// f = {3, 5, 4} it returns index 2; compile with
// -DTERNARY_LITERAL_REVERSE to see that reading fail.
template<class F>
int ternSearchRev(int a, int b, F f) {
	assert(a <= b);
	while (b - a >= 5) {
		int mid = a + (b - a) / 2;
		if (f(mid) <= f(mid+1)) a = mid; // (A)
		else b = mid+1;
	}
#ifdef TERNARY_LITERAL_REVERSE
	for (int i = b; i >= a+1; i--) if (f(a) < f(i)) a = i; // (B)
	return a;
#else
	for (int i = b-1; i >= a; i--) if (f(b) < f(i)) b = i; // (B)
	return b;
#endif
}
// Minimizing variant.
template<class F>
int ternSearchMin(int a, int b, F f) {
	assert(a <= b);
	while (b - a >= 5) {
		int mid = a + (b - a) / 2;
		if (f(mid) > f(mid+1)) a = mid; // (A)
		else b = mid+1;
	}
	rep(i,a+1,b+1) if (f(a) > f(i)) a = i; // (B)
	return a;
}

mt19937 rng(12345);
int rnd(int lo, int hi) { return uniform_int_distribution<int>(lo, hi)(rng); }

int main() {
	// random unimodal arrays: strictly increasing, then non-increasing
	rep(it,0,2000000) {
		int n = it < 1000000 ? rnd(1, 12) : rnd(1, 200);
		int off = rnd(-300, 300), peak = rnd(0, n-1);
		vi v(n);
		v[peak] = rnd(-5, 5);
		for (int i = peak-1; i >= 0; i--) v[i] = v[i+1] - rnd(1, 3);
		rep(i,peak+1,n) v[i] = v[i-1] - (rnd(0, 2) ? 0 : rnd(0, 3));
		int calls = 0;
		auto f = [&](int i) {
			assert(off <= i && i < off + n);
			calls++;
			return v[i - off];
		};
		int r = ternSearch(off, off + n - 1, f);
		assert(r == off + peak);
		assert(calls <= 2 * 12 + 12);
		// minimizing variant on negated values
		auto g = [&](int i) {
			assert(off <= i && i < off + n);
			return -v[i - off];
		};
		assert(ternSearchMin(off, off + n - 1, g) == off + peak);
		// reversed variant on the mirrored array: largest maximizer
		auto h = [&](int i) {
			assert(off <= i && i < off + n);
			return v[n - 1 - (i - off)];
		};
		assert(ternSearchRev(off, off + n - 1, h) == off + n - 1 - peak);
	}
	// large ranges, including ones where a + b does not fit in an int
	{
		const int B = 2000000000;
		vector<pii> ranges = {{0, B}, {-B, 0}, {B - 100, B}, {-B, -B + 100},
			{1000000000, B}, {-B, -1000000000}, {0, INT_MAX - 1},
			{INT_MIN, -1}, {INT_MAX - 7, INT_MAX - 1}, {-5, 5}, {7, 7}};
		for (auto [lo, hi] : ranges) rep(it,0,2000) {
			ll span = (ll)hi - lo;
			ll peak;
			if (it == 0) peak = lo;
			else if (it == 1) peak = hi;
			else peak = lo + (ll)(rng() % (unsigned ll)(span + 1));
			int calls = 0;
			auto f = [&](int i) {
				assert(lo <= i && i <= hi);
				calls++;
				return -abs((ll)i - peak);
			};
			assert(ternSearch(lo, hi, f) == peak);
			assert(calls <= 2 * 32 + 12);
		}
	}
	cout << "Tests passed!" << endl;
}
