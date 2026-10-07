#include "../utilities/template.h"

#include "../../content/strings/Manacher.h"

bool isPal(const string& s, int l, int r) { // [l, r)
	if (l < 0 || r > sz(s)) return false;
	while (l < r) if (s[l++] != s[--r]) return false;
	return true;
}

void test(const string& s) {
	int n = sz(s);
	array<vi, 2> p = manacher(s);
	assert(sz(p[0]) == n + 1 && sz(p[1]) == n);
	// p[0][i]: even palindromes s[i-k, i+k), centered just left of position i
	rep(i,0,n+1) {
		int k = 0;
		while (isPal(s, i - k - 1, i + k + 1)) k++;
		assert(p[0][i] == k);
	}
	// p[1][i]: odd palindromes s[i-k, i+k]
	rep(i,0,n) {
		int k = 0;
		while (isPal(s, i - k - 1, i + k + 2)) k++;
		assert(p[1][i] == k);
	}
}

template<class F>
void gen(string& s, int at, int alpha, F f) {
	if (at == sz(s)) f();
	else rep(i,0,alpha) {
		s[at] = (char)('a' + i);
		gen(s, at+1, alpha, f);
	}
}

int main() {
	test("");
	test("a");
	// all strings: binary up to length 16, ternary up to 10, 4 letters up to 8
	rep(n,0,17) { string s(n, 'x'); gen(s, 0, 2, [&]() { test(s); }); }
	rep(n,0,11) { string s(n, 'x'); gen(s, 0, 3, [&]() { test(s); }); }
	rep(n,0,9) { string s(n, 'x'); gen(s, 0, 4, [&]() { test(s); }); }
	// random medium strings, including bytes >= 128
	mt19937 rng(2024);
	rep(it,0,20000) {
		int n = (int)(rng() % 60), alpha = (int)(rng() % 3) + 1;
		string s(n, 'a');
		for (char& c : s) c = (char)(250 + rng() % alpha);
		test(s);
	}
	// structured: all equal, alternating, Fibonacci word
	rep(n,0,300) {
		test(string(n, 'a'));
		string s(n, 'a');
		rep(i,0,n) if (i & 1) s[i] = 'b';
		test(s);
	}
	string a = "a", b = "ab";
	while (sz(b) < 1500) { string c = b + a; a = b; b = c; }
	test(b);
	cout<<"Tests passed!"<<endl;
}
