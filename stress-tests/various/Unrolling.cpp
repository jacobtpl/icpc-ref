#include "../utilities/template.h"

// Unrolling.h is a code pattern, not a compilable header; this is
// the same pattern with a concrete body, checked against a plain loop.

mt19937 rng(31337);
int rnd(int lo, int hi) { return uniform_int_distribution<int>(lo, hi)(rng); }

const int OFF = 64;
int cnt[256], ord[256], no;

void unrolled(int from, int to) {
#define F {cnt[i + OFF]++; ord[no++] = i; ++i;}
	int i = from;
	while (i&3 && i < to) F // for alignment, if needed
	while (i + 4 <= to) { F F F F }
	while (i < to) F
#undef F
}

ll sumUnrolled(const int* a, int from, int to) {
	ll s = 0;
#define F {s += a[i]; ++i;}
	int i = from;
	while (i&3 && i < to) F
	while (i + 4 <= to) { F F F F }
	while (i < to) F
#undef F
	return s;
}

int main() {
	// exhaustive: every index in [from, to) visited once, in order
	rep(from,-OFF,OFF+1) rep(to,-OFF,OFF+1) {
		memset(cnt, 0, sizeof cnt); no = 0;
		unrolled(from, to);
		assert(no == max(0, to - from));
		rep(i,0,no) assert(ord[i] == from + i);
		rep(i,-OFF,OFF+1) assert(cnt[i + OFF] == (from <= i && i < to));
	}
	rep(it,0,200000) {
		int n = rnd(0, 100);
		vi a(n + 1);
		for (int& x : a) x = rnd(INT_MIN, INT_MAX);
		int from = rnd(0, n), to = rnd(0, n);
		ll s = 0;
		rep(i,from,to) s += a[i];
		assert(sumUnrolled(a.data(), from, to) == s);
	}
	cout << "Tests passed!" << endl;
}
