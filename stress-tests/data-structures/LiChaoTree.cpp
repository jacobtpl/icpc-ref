#include "../utilities/template.h"

#include "../../content/data-structures/LiChaoTree.h"
#include "../../content/data-structures/LineContainer.h"

typedef __int128 LL;

mt19937_64 rng(1234);
ll rnd(ll lo, ll hi) { // inclusive, hi - lo may exceed ll
	typedef unsigned long long ull;
	return (ll)((ull)lo + rng() % ((ull)hi - (ull)lo + 1));
}

struct Brute {
	ll lo, hi;
	vector<array<ll, 4>> v; // k, m, a, b
	Brute(ll lo, ll hi) : lo(lo), hi(hi) {}
	void add(ll k, ll m) { v.push_back({k, m, lo, hi}); }
	void add(ll k, ll m, ll a, ll b) { v.push_back({k, m, a, b}); }
	ll query(ll x) {
		LL res = LLONG_MAX;
		for (auto& [k, m, a, b] : v)
			if (a <= x && x < b) res = min(res, (LL)k * x + m);
		return (ll)res;
	}
};

// small ranges, every x checked after every op
void testSmall() {
	rep(it,0,20000) {
		ll len = rnd(1, 20), lo = rnd(-30, 30), hi = lo + len;
		LiChao t(lo, hi);
		Brute b(lo, hi);
		int mode = (int)rnd(0, 2); // 0 lines, 1 segments, 2 mixed
		int ops = (int)rnd(0, 12);
		for (ll x = lo; x < hi; x++) assert(t.query(x) == LLONG_MAX);
		rep(o,0,ops) {
			ll k = rnd(-10, 10), m = rnd(-50, 50);
			bool seg = mode == 1 || (mode == 2 && rnd(0, 1));
			if (!seg) t.add(k, m), b.add(k, m);
			else {
				ll a = rnd(lo - 5, hi + 5), bb = rnd(lo - 5, hi + 5);
				if (rnd(0, 3) == 0) bb = a + rnd(0, 1); // empty / single point
				t.add(k, m, a, bb), b.add(k, m, a, bb);
			}
			for (ll x = lo; x < hi; x++) assert(t.query(x) == b.query(x));
		}
	}
}

// n = 1 range
void testSingle() {
	rep(it,0,2000) {
		ll lo = rnd(-1000, 1000);
		LiChao t(lo, lo + 1);
		Brute b(lo, lo + 1);
		assert(t.query(lo) == LLONG_MAX);
		rep(o,0,5) {
			ll k = rnd(-100, 100), m = rnd(-1000, 1000);
			if (rnd(0, 1)) t.add(k, m), b.add(k, m);
			else {
				ll a = rnd(lo - 2, lo + 2), bb = rnd(lo - 2, lo + 2);
				t.add(k, m, a, bb), b.add(k, m, a, bb);
			}
			assert(t.query(lo) == b.query(lo));
		}
	}
}

// big ranges and values near the ll limits
void testBig() {
	// |x| <= 2e9, |k| <= 2e9, |m| <= 5e18 -> |kx+m| <= 9e18
	rep(it,0,200) {
		ll lo = -2000000000, hi = 2000000001;
		LiChao t(lo, hi);
		Brute b(lo, hi);
		vector<ll> xs = {lo, lo + 1, hi - 2, hi - 1, 0, -1, 1};
		rep(i,0,20) xs.push_back(rnd(lo, hi - 1));
		rep(o,0,30) {
			ll k = rnd(-2000000000, 2000000000);
			ll m = rnd(-5000000000000000000LL, 5000000000000000000LL);
			if (rnd(0, 1)) t.add(k, m), b.add(k, m);
			else {
				ll a = rnd(lo, hi), bb = rnd(lo, hi);
				t.add(k, m, a, bb), b.add(k, m, a, bb);
				xs.push_back(a); xs.push_back(bb - 1);
			}
			for (ll x : xs) if (lo <= x && x < hi)
				assert(t.query(x) == b.query(x));
		}
	}
	// range of width ~8e18, slopes in -1..1
	rep(it,0,100) {
		ll lo = it % 2 ? -4000000000000000000LL : 0;
		ll hi = 4000000000000000000LL;
		LiChao t(lo, hi);
		Brute b(lo, hi);
		vector<ll> xs = {lo, lo + 1, hi - 2, hi - 1, (lo + hi) / 2};
		rep(i,0,20) xs.push_back(rnd(lo, hi - 1));
		rep(o,0,30) {
			ll k = rnd(-1, 1), m = rnd(-1000000000000000000LL, 1000000000000000000LL);
			if (rnd(0, 1)) t.add(k, m), b.add(k, m);
			else {
				ll a = rnd(lo, hi), bb = rnd(lo, hi);
				t.add(k, m, a, bb), b.add(k, m, a, bb);
				xs.push_back(a); xs.push_back(bb - 1);
			}
			for (ll x : xs) if (lo <= x && x < hi)
				assert(t.query(x) == b.query(x));
		}
	}
	// constants at the limits
	{
		LiChao t(-5, 5);
		t.add(0, LLONG_MAX);
		rep(x,-5,5) assert(t.query(x) == LLONG_MAX);
		t.add(0, LLONG_MIN, 0, 3);
		rep(x,-5,5) assert(t.query(x) == (0 <= x && x < 3 ? LLONG_MIN : LLONG_MAX));
		t.add(1, 0);
		rep(x,-5,5) assert(t.query(x) == (0 <= x && x < 3 ? LLONG_MIN : x));
	}
}

// max variant vs LineContainer
void testMax() {
	rep(it,0,300) {
		ll lo = -1000000000, hi = 1000000001;
		LiChao t(lo, hi);
		LineContainer lc;
		rep(o,0,200) {
			ll k = rnd(-1000000000, 1000000000);
			ll m = rnd(-1000000000000000000LL, 1000000000000000000LL);
			t.add(-k, -m);
			lc.add(k, m);
			rep(q,0,5) {
				ll x = rnd(lo, hi - 1);
				assert(-t.query(x) == lc.query(x));
			}
		}
	}
}

// duplicate / parallel lines on a mid-size range
void testParallel() {
	ll lo = -1000, hi = 1000;
	LiChao t(lo, hi);
	Brute b(lo, hi);
	rep(o,0,300) {
		ll k = rnd(-2, 2), m = rnd(-3, 3);
		if (o % 3 == 0) t.add(k, m), b.add(k, m);
		else {
			ll a = rnd(lo, hi), bb = rnd(lo, hi);
			t.add(k, m, a, bb), b.add(k, m, a, bb);
		}
	}
	for (ll x = lo; x < hi; x++) assert(t.query(x) == b.query(x));
}

// speed: 2e5 random lines + segments on C = 2^31, 2e5 queries
void testSpeed() {
	ll lo = -(1LL << 30), hi = 1LL << 30;
	LiChao t(lo, hi);
	ll sink = 0;
	rep(i,0,200000) {
		ll k = rnd(-1000000, 1000000), m = rnd(-1000000000000LL, 1000000000000LL);
		if (i % 2) t.add(k, m);
		else {
			ll a = rnd(lo, hi), bb = rnd(lo, hi);
			t.add(k, m, a, bb);
		}
	}
	rep(i,0,200000) sink ^= t.query(rnd(lo, hi - 1));
	assert(sink != 12345); // keep the queries alive
}

int main() {
	testSmall();
	testSingle();
	testBig();
	testMax();
	testParallel();
	testSpeed();
	cout << "Tests passed!" << endl;
}
