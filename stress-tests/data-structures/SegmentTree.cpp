#include "../utilities/template.h"

static unsigned RA = 1231231;
int ra() {
	RA *= 574841;
	RA += 14;
	return RA >> 1;
}

// The header can only be included once (#pragma once), so to also run the real
// header with a non-commutative f, its max() is routed through segF, which is
// plain max unless segMode is set (see nonabelianF below).
int segMode = 0;
int segF(int a, int b);

namespace maximum {

#define max(a, b) segF(a, b)
#include "../../content/data-structures/SegmentTree.h"
#undef max

}

namespace nonabelian {

// https://en.wikipedia.org/wiki/Dihedral_group_of_order_6
const int lut[6][6] = {
	{0, 1, 2, 3, 4, 5},
	{1, 0, 4, 5, 2, 3},
	{2, 5, 0, 4, 3, 1},
	{3, 4, 5, 0, 1, 2},
	{4, 3, 1, 2, 5, 0},
	{5, 2, 3, 1, 0, 4}
};

struct Tree {
	typedef int T;
	const T unit = 0;
	T f(T a, T b) { return lut[a][b]; }
	vector<T> s; int n;
	Tree(int n = 0, T def = 0) : s(2*n, def), n(n) {}
	void update(int pos, T val) {
		for (s[pos += n] = val; pos > 1; pos /= 2)
			s[pos / 2] = f(s[pos & ~1], s[pos | 1]);
	}
	T query(int b, int e) { // query [b, e)
		T ra = unit, rb = unit;
		for (b += n, e += n; b < e; b /= 2, e /= 2) {
			if (b % 2) ra = f(ra, s[b++]);
			if (e % 2) rb = f(s[--e], rb);
		}
		return f(ra, rb);
	}
};

}

// segMode = 1: a non-commutative monoid on {INT_MIN} + D6 in which INT_MIN
// (the header's unit) is the identity.
int nonabelianF(int a, int b) {
	return a == INT_MIN ? b : b == INT_MIN ? a : nonabelian::lut[a][b];
}
int segF(int a, int b) { return segMode ? nonabelianF(a, b) : max(a, b); }

void testMore() {
	mt19937 rng(7);
	auto ri = [&](int a, int b) { return uniform_int_distribution<int>(a, b)(rng); };
	// max-tree: all sizes up to 70 (non powers of two), negative and extreme values
	rep(n,0,71) rep(rounds,0,20) {
		int mode = ri(0, 2);
		auto val = [&]() {
			return mode == 0 ? ri(INT_MIN, INT_MAX) : mode == 1 ? ri(-3, 3) :
				(ri(0, 1) ? INT_MIN : INT_MAX);
		};
		int def = ri(0, 1) ? maximum::Tree::unit : val();
		maximum::Tree tr(n, maximum::Tree::unit);
		vi v(n, maximum::Tree::unit);
		if (def != maximum::Tree::unit) // a non-unit default needs explicit updates
			rep(i,0,n) tr.update(i, def), v[i] = def;
		rep(it,0,300) {
			if (n && ri(0, 2)) {
				int i = ri(0, n - 1);
				tr.update(i, v[i] = val());
			}
			int i = ri(0, n), j = ri(0, n);
			int ma = INT_MIN;
			rep(k,i,j) ma = max(ma, v[k]);
			assert(tr.query(i, j) == ma); // i >= j gives unit
		}
	}
	// non-commutative f on the real header
	segMode = 1;
	rep(n,1,40) rep(rounds,0,20) {
		maximum::Tree tr(n);
		vi v(n, INT_MIN);
		rep(it,0,300) {
			if (ri(0, 2)) {
				int i = ri(0, n - 1);
				tr.update(i, v[i] = ri(0, 5));
			}
			int i = ri(0, n), j = ri(0, n);
			int r = INT_MIN;
			rep(k,i,j) r = nonabelianF(r, v[k]);
			assert(tr.query(i, j) == r);
		}
	}
	segMode = 0;
}

void bench() {
	for (int N : {200000, 1000000}) {
		maximum::Tree tr(N);
		auto t0 = chrono::steady_clock::now();
		ll sum = 0;
		rep(i,0,N) tr.update(i, ra());
		rep(it,0,2000000) {
			tr.update(ra() % N, ra());
			int i = ra() % N, j = ra() % N;
			if (i > j) swap(i, j);
			sum += tr.query(i, j+1);
		}
		cerr << "N=" << N << ", 2e6 updates + 2e6 queries: " << chrono::duration<double>(
			chrono::steady_clock::now() - t0).count() << " s (" << sum << ")\n";
	}
}

int main(int argc, char**) {
	if (argc > 1) return bench(), 0; // benchmark mode: ./a.out bench
	testMore();

	{
		maximum::Tree t(0);
		assert(t.query(0, 0) == t.unit);
	}

	if (0) {
		const int N = 10000;
		maximum::Tree tr(N);
		ll sum = 0;
		rep(it,0,1000000) {
			tr.update(ra() % N, ra());
			int i = ra() % N;
			int j = ra() % N;
			if (i > j) swap(i, j);
			int v = tr.query(i, j+1);
			sum += v;
		}
		cout << sum << endl;
		// return 0;
	}

	rep(n,1,10) {
		maximum::Tree tr(n);
		vi v(n, maximum::Tree::unit);
		rep(it,0,1000000) {
			int i = rand() % (n+1), j = rand() % (n+1);
			int x = rand() % (n+2);

			int r = rand() % 100;
			if (r < 30) {
				int ma = tr.unit;
				rep(k,i,j) ma = max(ma, v[k]);
				assert(ma == tr.query(i,j));
			}
			else {
				i = min(i, n-1);
				tr.update(i, x);
				v[i] = x;
			}
		}
	}

	rep(n,1,10) {
		nonabelian::Tree tr(n);
		vi v(n);
		rep(it,0,1000000) {
			int i = rand() % (n+1), j = rand() % (n+1);
			int x = rand() % 6;

			int r = rand() % 100;
			if (r < 30) {
				int ma = tr.unit;
				rep(k,i,j) ma = nonabelian::lut[ma][v[k]];
				assert(ma == tr.query(i,j));
			}
			else {
				i = min(i, n-1);
				tr.update(i, x);
				v[i] = x;
			}
		}
	}
	cout<<"Tests passed!"<<endl;
}
