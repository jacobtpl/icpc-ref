#include "../utilities/template.h"

#include "../../content/strings/KMP.h"

template<class F>
void gen(string& s, int at, int alpha, F f) {
	if (at == sz(s)) f();
	else {
		rep(i,0,alpha) {
			s[at] = (char)('a' + i);
			gen(s, at+1, alpha, f);
		}
	}
}

void test(const string& s) {
	vi p = pi(s);
	rep(i,0,sz(s)) {
		int maxlen = -1;
		rep(len,0,i+1) {
			rep(j,0,len) {
				if (s[j] != s[i+1 - len + j]) goto fail;
			}
			maxlen = len;
fail:;
		}
		assert(maxlen == p[i]);
	}
}

// Occurrences of pat in s, by definition.
vi brute(const string& s, const string& pat) {
	vi res;
	rep(i,0,sz(s)-sz(pat)+1)
		if (s.compare(i, sz(pat), pat) == 0) res.push_back(i);
	return res;
}

void testMatch(const string& s, const string& pat) {
	vi want = brute(s, pat);
	assert(match(s, pat) == want);
	assert(match2(s, pat) == want);
}

mt19937 rng(12345);
string randStr(int n, int alpha) {
	string s(n, 'a');
	for (char& c : s) c = (char)('a' + rng() % alpha);
	return s;
}

void testMatches() {
	// edge cases: single characters, pattern == text, pattern longer than text
	testMatch("a", "a");
	testMatch("a", "b");
	testMatch("", "a");
	testMatch("ab", "abc");
	testMatch("aaaa", "a");
	testMatch("aaaa", "aa");
	testMatch("abab", "abab");
	testMatch("abcab", "ab");
	// exhaustive: every text of length <= 7 and pattern of length 1..4 over {a,b}
	rep(n,0,8) rep(m,1,5) rep(a,0,1<<n) rep(b,0,1<<m) {
		string s(n, 'a'), pat(m, 'a');
		rep(i,0,n) if (a >> i & 1) s[i] = 'b';
		rep(i,0,m) if (b >> i & 1) pat[i] = 'b';
		testMatch(s, pat);
	}
	// random, small alphabets, patterns frequently taken from the text
	rep(it,0,200000) {
		int alpha = (int)(rng() % 3) + 1, n = (int)(rng() % 25);
		string s = randStr(n, alpha), pat;
		if (n && rng() % 2) {
			int a = (int)(rng() % n), len = (int)(rng() % (n - a)) + 1;
			pat = s.substr(a, len);
		} else pat = randStr((int)(rng() % 6) + 1, alpha);
		testMatch(s, pat);
	}
	// full byte range except NUL, which match() uses as a separator
	rep(it,0,20000) {
		int n = (int)(rng() % 20);
		string s(n, 'a');
		for (char& c : s) c = (char)(rng() % 2 ? 200 + rng() % 3 : 1 + rng() % 2);
		int m = (int)(rng() % 3) + 1;
		string pat = n >= m && rng() % 2 ? s.substr(rng() % (n - m + 1), m) : string(m, (char)200);
		testMatch(s, pat);
	}
	// large periodic inputs
	testMatch(string(200000, 'a'), string(1000, 'a'));
	testMatch(string(200000, 'a'), string(999, 'a') + "b");
	{
		string s;
		rep(i,0,50000) s += "ab";
		testMatch(s, "abab");
		testMatch(s, "ba");
	}
#ifdef TEST_EMPTY_PATTERN
	// Undocumented precondition: an empty pattern is not supported
	// (both functions return 1..n instead of 0..n).
	testMatch("abc", "");
#endif
}

int main() {
	testMatches();
	// string str; cin >> str; for(auto &x: pi(str)) cout << x; cout << endl;
	// test ~3^12 strings
	rep(n,0,13) {
		string s(n, 'x');
		gen(s, 0, 3, [&]() {
			test(s);
		});
	}
	// then ~4^10 strings
	rep(n,0,11) {
		string s(n, 'x');
		gen(s, 0, 4, [&]() {
			test(s);
		});
	}
	cout<<"Tests passed!"<<endl;
}
