#include "../utilities/template.h"

#include "../../content/strings/Lyndon.h"
#include "../../content/strings/MinRotation.h"

template<class S> bool isLyndon(const S& w) {
	int n = sz(w);
	if (n == 0) return false;
	rep(i,1,n) {
		S suf(w.begin() + i, w.end());
		if (!(w < suf)) return false;
	}
	return true;
}

// Repeatedly take the longest Lyndon prefix.
template<class S> vi lyndonSlow(const S& s) {
	int n = sz(s), i = 0;
	vi r = {0};
	while (i < n) {
		int best = i + 1;
		rep(j,i+1,n+1) if (isLyndon(S(s.begin() + i, s.begin() + j))) best = j;
		r.push_back(i = best);
	}
	return r;
}

template<class S> void check(const S& s) {
	int n = sz(s);
	vi r = lyndon(s);
	assert(r == lyndonSlow(s));
	assert(sz(r) >= 1 && r[0] == 0 && r.back() == n);
	rep(i,1,sz(r)) {
		S w(s.begin() + r[i-1], s.begin() + r[i]);
		assert(isLyndon(w));
		if (i >= 2) assert(!(S(s.begin() + r[i-2], s.begin() + r[i-1]) < w));
	}
	if (n > 0) {
		int st = r[sz(r) - 2], bst = 0;
		rep(i,1,n) if (S(s.begin() + i, s.end()) < S(s.begin() + bst, s.end())) bst = i;
		assert(st == bst);
		S ss(s); ss.insert(ss.end(), all(s));
		vi r2 = lyndon(ss);
		int rot = -1;
		for (int x : r2) if (x < n) rot = x;
		S best(ss.begin() + rot, ss.begin() + rot + n);
		rep(i,0,n) assert(!(S(ss.begin() + i, ss.begin() + i + n) < best));
	}
}

void testPerf() {
	string s;
	rep(i,0,10000000) s += (char)('a' + (rand() % 400000 < 2));
	vi r = lyndon(s);
	cout << sz(r) << endl;
	s.assign(10000000, 'a');
	r = lyndon(s);
	assert(sz(r) == 10000001);
}

int main() {
	check(string());
	check(vi());
	check(string("a"));
	check(string("abaab"));
	assert(lyndon(string("abaab")) == (vi{0, 2, 5}));
	assert(lyndon(string()) == (vi{0}));
	assert(lyndon(vi{5}) == (vi{0, 1}));
	assert(lyndon(vi{3, 2, 1}) == (vi{0, 1, 2, 3}));
	assert(lyndon(vi{1, 2, 3}) == (vi{0, 3}));
	assert(lyndon(vi{2, 2, 2}) == (vi{0, 1, 2, 3}));
	rep(it,0,200000) {
		int n = rand() % 12, a = rand() % 4 + 1;
		string s;
		rep(i,0,n) s += (char)('a' + rand() % a);
		check(s);
		if (n > 0) {
			string ss = s + s;
			vi r = lyndon(ss);
			int rot = -1;
			for (int x : r) if (x < n) rot = x;
			rotate(s.begin(), s.begin() + rot, s.end());
			rotate(ss.begin(), ss.begin() + minRotation(ss) % n, ss.end());
			assert(s == ss.substr(0, n));
		}
	}
	rep(it,0,100000) {
		int n = rand() % 10;
		vi v;
		rep(i,0,n) {
			int t = rand() % 6;
			v.push_back(t == 0 ? INT_MIN : t == 1 ? INT_MAX : rand() % 3 - 1);
		}
		check(v);
		vector<ll> w(all(v));
		for (ll& x : w) x = x == INT_MIN ? LLONG_MIN : x == INT_MAX ? LLONG_MAX : x;
		check(w);
		vector<signed char> c;
		for (int x : v) c.push_back((signed char)(x == INT_MIN ? -128 : x == INT_MAX ? 127 : x));
		check(c);
	}
	rep(it,0,2000) {
		int n = rand() % 60 + 1, p = rand() % 5 + 1;
		string s;
		rep(i,0,n) s += (char)('a' + (i % p == 0));
		if (rand() % 2) s[rand() % n] = (char)('a' + rand() % 3);
		check(s);
	}
	testPerf();
	cout << "Tests passed!" << endl;
}
