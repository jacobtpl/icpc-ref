#include "../utilities/template.h"
// from content/contest/template.cpp
#define pb push_back
template<class A, class B>
bool ckmin(A &a, const B& b) { return b < a ? a = b, 1 : 0; }

#include "../../content/strings/PalTree.h"

bool isPal(const string& s, int l, int r) { // [l, r)
	for (r--; l < r; l++, r--) if (s[l] != s[r]) return 0;
	return 1;
}

void test(const string& s, bool online) {
	const int inf = PalTree::INF;
	int n = sz(s);
	PalTree pt;
	set<string> pals;
	// f[i][b] = min #palindromes with parity b partitioning s[0, i)
	vector<array<int,2>> f(n + 1, {inf, inf});
	f[0][0] = 0;
	rep(i,1,n+1) {
		pt.addChar(s[i-1]);
		int longest = 0;
		rep(j,0,i) if (isPal(s, j, i)) {
			pals.insert(s.substr(j, i - j));
			longest = max(longest, i - j);
			rep(b,0,2) if (f[j][b^1] < inf)
				f[i][b] = min(f[i][b], f[j][b^1] + 1);
		}
		if (online || i == n) {
			assert(pt.d[pt.last].len == longest);
			assert(sz(pt.d) - 2 == sz(pals));
			assert(sz(pt.ans) == i + 1);
			rep(k,0,i+1) assert(pt.ans[k] == f[k]);
		}
	}
	assert(pt.ans == f);
	assert(pt.d[0].len == 0 && pt.d[1].len == -1);
	// recover the palindrome of every node through the trie
	vector<string> str(sz(pt.d));
	vi seen(sz(pt.d));
	seen[0] = seen[1] = 1;
	rep(v,0,sz(pt.d)) {
		assert(seen[v]); // children are created after their parent
		rep(c,0,PalTree::ASZ) if (int u = pt.d[v].to[c]) {
			assert(u > v && u >= 2 && !seen[u]);
			seen[u] = 1;
			string ch(1, (char)('a' + c));
			str[u] = v == 1 ? ch : ch + str[v] + ch;
			assert(pt.d[u].len == sz(str[u]));
		}
	}
	set<string> got(str.begin() + 2, str.end());
	assert(got == pals);
	map<string,int> id;
	rep(v,2,sz(pt.d)) id[str[v]] = v;
	id[""] = 0;
	pt.resolveOc();
	rep(v,2,sz(pt.d)) {
		const string& p = str[v];
		int m = sz(p), oc = 0, lk = 0;
		rep(i,0,n-m+1) oc += s.compare(i, m, p) == 0;
		assert(pt.d[v].oc == oc);
		rep(k,1,m) if (isPal(p, k, m)) { lk = id.at(p.substr(k)); break; }
		assert(pt.d[v].link == lk);
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
	test("", 1);
	rep(c,0,26) test(string(1, (char)('a' + c)), 1);
	// exhaustive small
	rep(n,0,13) { string s(n, 'x'); gen(s, 0, 2, [&]() { test(s, 1); }); }
	rep(n,0,9) { string s(n, 'x'); gen(s, 0, 3, [&]() { test(s, 1); }); }
	rep(n,0,7) { string s(n, 'x'); gen(s, 0, 4, [&]() { test(s, 1); }); }
	mt19937 rng(2024);
	auto rnd = [&](int a, int b) {
		return (int)(rng() % (unsigned)(b - a + 1)) + a; };
	rep(it,0,6000) {
		int n = rnd(1, it % 20 == 0 ? 300 : 40), kind = rnd(0, 4);
		int alpha = kind == 0 ? 26 : rnd(1, 3);
		string s(n, 'a');
		rep(i,0,n) s[i] = (char)('a' + rnd(0, alpha-1));
		if (kind == 2) { int p = rnd(1, 6); rep(i,p,n) s[i] = s[i-p]; }
		if (kind == 3) rep(i,0,n) s[i] = (char)('a' + __builtin_ctz(i+1));
		if (kind == 4) rep(i,0,n) s[i] = (char)('y' + (rnd(0, 15) == 0));
		test(s, n <= 40);
	}
	// all-equal and fibonacci-like strings at larger sizes
	test(string(400, 'z'), 0);
	string a = "a", b = "ab";
	while (sz(b) < 400) { a = b + a; swap(a, b); }
	test(b, 0);
	cout<<"Tests passed!"<<endl;
}
