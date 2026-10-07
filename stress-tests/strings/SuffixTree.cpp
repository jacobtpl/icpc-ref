#include "../utilities/template.h"

#include "../../content/strings/SuffixTree.h"

// The struct is ~24 MB, so it must not live on the stack.
SuffixTree* build(const string& a) { return new SuffixTree(a); }

set<string> substrings(const string& s) {
	set<string> r;
	rep(i,0,sz(s)) rep(l,1,sz(s)-i+1) r.insert(s.substr(i, l));
	return r;
}

int leaves, inner;
void dfs(SuffixTree& st, const string& a, int v, string& cur,
		set<string>& out, bool complete) {
	int n = sz(a), kids = 0;
	rep(c,0,SuffixTree::ALPHA) if (st.t[v][c] != -1) {
		int u = st.t[v][c];
		kids++;
		assert(2 <= u && u < st.m);
		assert(st.p[u] == v);
		assert(0 <= st.l[u] && st.l[u] < st.r[u] && st.r[u] <= n);
		assert(st.toi(a[st.l[u]]) == c);
		size_t before = cur.size();
		rep(i,st.l[u],st.r[u]) cur += a[i], out.insert(cur);
		dfs(st, a, u, cur, out, complete);
		cur.resize(before);
	}
	if (!kids) {
		assert(st.r[v] == n);
		leaves++;
		if (complete) // leaf = suffix
			assert(a.compare(n - sz(cur), sz(cur), cur) == 0);
	} else if (v) {
		inner++;
		assert(kids >= 2);
		// suffix link drops the first character
		int u = st.s[v], len = sz(cur) - 1;
		string w;
		while (u > 0) {
			w = a.substr(st.l[u], st.r[u] - st.l[u]) + w;
			u = st.p[u];
		}
		assert(u == 0 && w == cur.substr(1) && sz(w) == len);
	}
}

// complete: the last character of a is unique
void test(const string& a, bool complete) {
	SuffixTree* st = build(a);
	assert(st->l[0] == -1 && st->r[0] == 0);
	assert(st->m <= 2 * sz(a) + 2);
	set<string> got;
	string cur;
	leaves = inner = 0;
	if (!a.empty()) dfs(*st, a, 0, cur, got, complete);
	assert(got == substrings(a));
	if (complete) assert(leaves == sz(a));
	delete st;
}

// lcs() as used by LCS(), with separators inside a..z so that the
// documented ALPHA = 28 change is not needed
void testLcs(const string& s, const string& t) {
	string a = s + 'y' + t + 'z';
	SuffixTree* st = build(a);
	st->best = {0, 0};
	st->lcs(0, sz(s), sz(s) + 1 + sz(t), 0);
	int len = st->best.first, pos = st->best.second, exp = 0;
	rep(i,0,sz(s)) rep(j,0,sz(t)) {
		int k = 0;
		while (i + k < sz(s) && j + k < sz(t) && s[i+k] == t[j+k]) k++;
		exp = max(exp, k);
	}
	assert(len == exp);
	if (len) {
		assert(0 <= pos && pos + len <= sz(a));
		string w = a.substr(pos, len);
		assert(s.find(w) != string::npos && t.find(w) != string::npos);
	}
	delete st;
}

int main() {
#ifdef LCS_STATIC // opt-in: LCS() puts a ~24 MB SuffixTree on the stack
	cout << SuffixTree::LCS("ab", "b").first << endl;
#endif
	mt19937 rng(99);
	auto rnd = [&](int a, int b) {
		return (int)(rng() % (unsigned)(b - a + 1)) + a; };
	test("", 0); test("a", 0); test("z", 1); test("aa", 0);
	rep(it,0,700) {
		int n = rnd(0, it % 10 == 0 ? 60 : 14), kind = rnd(0, 3);
		int alpha = kind == 0 ? 25 : rnd(1, 3);
		string s(n, 'a');
		rep(i,0,n) s[i] = (char)('a' + rnd(0, alpha-1));
		if (kind == 2) { int p = rnd(1, 4); rep(i,p,n) s[i] = s[i-p]; }
		if (it % 2) test(s + 'z', 1);
		else test(s, 0);
	}
	rep(it,0,500) {
		int alpha = rnd(1, 3), n = rnd(0, 12), m = rnd(0, 12);
		string s, t;
		rep(i,0,n) s += (char)('a' + rnd(0, alpha-1));
		rep(i,0,m) t += (char)('a' + rnd(0, alpha-1));
		if (n && rnd(0, 2) == 0) t += s.substr(rnd(0, n-1));
		testLcs(s, t);
	}
	{ // documented maximum length (N ~ 2*maxlen+10)
		string s(100000, 'a');
		rep(i,0,sz(s)-1) s[i] = (char)('a' + rnd(0, 1));
		s.back() = 'z';
		SuffixTree* st = build(s);
		assert(st->m <= SuffixTree::N);
		int lv = 0;
		rep(v,2,st->m) {
			bool leaf = 1;
			rep(c,0,26) if (st->t[v][c] != -1) leaf = 0;
			lv += leaf;
		}
		assert(lv == sz(s));
		delete st;
	}
	cout<<"Tests passed!"<<endl;
}
