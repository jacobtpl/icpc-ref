#include "../utilities/template.h"

#include "../../content/combinatorial/multinomial.h"

// exact oracle via Legendre's formula; returns -1 if result > LLONG_MAX
ll oracle(const vi& v) {
	int s = accumulate(all(v), 0);
	auto leg = [](int n, int p) { int e = 0; while (n) e += n /= p; return e; };
	unsigned __int128 r = 1;
	rep(p,2,s+1) {
		bool pr = 1;
		for (int d = 2; d * d <= p; d++) if (p % d == 0) pr = 0;
		if (!pr) continue;
		int e = leg(s, p);
		for (int x : v) e -= leg(x, p);
		assert(e >= 0);
		while (e--) { r *= p; if (r > (unsigned __int128)LLONG_MAX) return -1; }
	}
	return (ll)r;
}

int tested = 0;
void check(vi v) {
	ll want = oracle(v);
	if (want < 0) return; // result does not fit in ll
	tested++;
	ll got = multinomial(v);
	if (got != want) {
		cout << "multinomial(";
		for (int x : v) cout << x << ' ';
		cout << ") = " << got << ", expected " << want << endl;
		exit(1);
	}
}

int main() {
	mt19937 rng(2024);
	auto rnd = [&](int a, int b) { return uniform_int_distribution<int>(a, b)(rng); };
	check({}); check({0}); check({7}); check({0,0,0}); check({1,1,1,1});
	check({100000}); check({100000, 1}); check({1, 100000});
	check({100000, 3}); check({3, 100000}); check({0, 100000, 0, 2});
	// 20! is the largest factorial fitting in ll
	check(vi(20, 1));
	// exhaustive small
	rep(a,0,13) rep(b,0,13) rep(c,0,13) { check({a,b,c}); check({a,b}); }
	rep(a,0,7) rep(b,0,7) rep(c,0,7) rep(d,0,7) rep(e,0,7) check({a,b,c,d,e});
	// binomials: all C(a+b, a) with a+b <= 66 (largest row where all fit is 66)
	rep(a,0,67) rep(b,0,67) check({a,b});
	// results close to LLONG_MAX
	check({31,31}); check({33,33}); check({30,36}); check({1,31,31});
	rep(it,0,300000) {
		int n = rnd(1, 8), s = rnd(0, 70);
		vi v(n);
		rep(i,0,s) v[rnd(0, n-1)]++;
		if (rnd(0, 3) == 0) sort(all(v));
		if (rnd(0, 3) == 0) sort(all(v), greater<int>());
		check(v);
	}
	assert(tested > 100000);
	cout << "Tests passed!" << endl;
}
