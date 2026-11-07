#include "../utilities/template.h"

#include "../../content/strings/SuffixAutomaton.h"

// Ordering used by map<char> (plain char), which is what kth() sorts by.
bool charLess(const string& a, const string& b) {
	return lexicographical_compare(all(a), all(b));
}

int walk(SA& sa, const string& s) {
	int p = 0;
	for (char c : s) {
		auto it = sa.adj[p].find(c);
		if (it == sa.adj[p].end()) return -1;
		p = it->second;
	}
	return p;
}

int countOcc(const string& s, const string& w) {
	if (w.empty()) return sz(s) + 1;
	int r = 0;
	for (size_t i = s.find(w); i != string::npos; i = s.find(w, i + 1)) r++;
	return r;
}

string randString(int n, int alpha, bool bytes) {
	string s;
	rep(i,0,n) s += bytes ? (char)(rand() % 256) : (char)('a' + rand() % alpha);
	return s;
}

void check(const vector<string>& strs, const vector<string>& probes) {
	SA sa;
	vector<int> ends;
	set<string> subs;
	for (auto& s : strs) {
		int p = 0;
		for (char c : s) p = sa.append(p, c);
		ends.push_back(p);
		rep(i,0,sz(s)) rep(j,i+1,sz(s)+1) subs.insert(s.substr(i, j - i));
	}
	int N = sa.N;
	assert(sa.link[0] == -1 && sa.dis[0] == 0);
	rep(i,1,N) assert(sa.link[i] >= 0 && sa.dis[sa.link[i]] < sa.dis[i]);
	rep(i,0,sz(strs)) assert(walk(sa, strs[i]) == ends[i]);

	// order: permutation sorted by dis
	vi o = sa.order();
	assert(sz(o) == N);
	vi seen(N);
	rep(i,0,N) assert(0 <= o[i] && o[i] < N && !seen[o[i]]++);
	rep(i,1,N) assert(sa.dis[o[i-1]] <= sa.dis[o[i]]);

	// longest string of every state, by processing states in order of dis
	vector<string> lng(N);
	vector<bool> got(N); got[0] = 1;
	for (int u : o) {
		assert(got[u]);
		for (auto [c, v] : sa.adj[u]) {
			assert(sa.dis[v] > sa.dis[u]);
			if (sa.dis[v] == sa.dis[u] + 1) lng[v] = lng[u] + c, got[v] = 1;
		}
	}
	// every state's strings are exactly the suffixes of lng with length in
	// (dis[link], dis]; they are distinct substrings and cover everything
	ll total = 0;
	rep(v,1,N) {
		int lo = sa.dis[sa.link[v]], hi = sa.dis[v];
		assert(sz(lng[v]) == hi);
		for (int len = lo + 1; len <= hi; len++) {
			string w = lng[v].substr(hi - len);
			assert(walk(sa, w) == v);
			assert(subs.count(w));
			total++;
		}
		string shorter = lng[v].substr(hi - lo);
		assert(walk(sa, shorter) == sa.link[v]);
	}
	assert(total == sz(subs));
	assert(sa.distinct() == sz(subs));

	// membership
	for (auto& w : subs) assert(walk(sa, w) != -1);
	for (auto& w : probes) if (!w.empty())
		assert((walk(sa, w) != -1) == (bool)subs.count(w));

	// kth vs sorted (by plain char) list of substrings
	vector<string> sorted(all(subs));
	sort(all(sorted), charLess);
	rep(k,0,sz(sorted)) assert(sa.kth(k) == sorted[k]);

	// occ vs naive counting, for each added string
	for (auto& s : strs) {
		vi cnt = sa.occ(s);
		assert(sz(cnt) == N);
		rep(v,1,N) {
			int lo = sa.dis[sa.link[v]], hi = sa.dis[v];
			for (int len = lo + 1; len <= hi; len++)
				assert(cnt[v] == countOcc(s, lng[v].substr(hi - len)));
		}
	}

	// lcs vs brute force
	for (auto& t : probes) {
		pii r = sa.lcs(t);
		int best = 0;
		rep(i,0,sz(t)) rep(j,i+1,sz(t)+1)
			if (subs.count(t.substr(i, j - i))) best = max(best, j - i);
		assert(r.first == best);
		assert(0 <= r.second && r.second + r.first <= sz(t));
		if (best) assert(subs.count(t.substr(r.second, r.first)));
	}
}

void testEdge() {
	// nothing added
	{
		SA sa;
		assert(sa.N == 1 && sa.distinct() == 0);
		assert(sa.order() == vi{0});
		assert(sa.lcs("abc") == pii(0, 3));
		assert(sa.lcs("") == pii(0, 0));
		assert(sa.occ("") == vi{0});
	}
	check({""}, {"", "a"});
	check({"", "", ""}, {"", "a"});
	check({"a"}, {"", "a", "b", "aa"});
	check({"a", "a", "a"}, {"", "a", "b", "aa"});
	check({"ab", "ab"}, {"ab", "ba", "b", "abab"});
	check({string(12, 'a')}, {string(13, 'a'), "aab", "baa"});
	check({string(6, 'a') + "b" + string(6, 'a')}, {"aba", "bab", "aaaaaaa"});
	check({"abcabcabc"}, {"cab", "bca", "abcabcabca"});
	check({"ab", "b", "a", "ba", "abab"}, {"aba", "bab", "ba", "c"});
	check({"abcdef", "cde", "abc", "def"}, {"abcdef", "cdef", "xabc"});
	// bytes >= 128 (negative char) and 0 bytes
	check({string("\x80\xff\x00\x7f\x80", 5)}, {string("\xff\x00", 2), "\x7f", "x"});
	check({string("\x00\x00\x01", 3), string("\x01\x00", 2)}, {string("\x00\x01\x00", 3), string("\x00", 1)});
}

void testOverflow() {
	// distinct() exceeds int, and kth() works near the top of the range.
	srand(7);
	int n = 70000;
	string s = randString(n, 3, false);
	SA sa;
	int p = 0;
	for (char c : s) p = sa.append(p, c);
	assert(sa.N <= 2 * n && sa.N > n);
	vi suf(n); iota(all(suf), 0);
	sort(all(suf), [&](int a, int b) {
		return string_view(s).substr(a) < string_view(s).substr(b);
	});
	ll distinct = (ll)n * (n + 1) / 2;
	rep(i,1,n) {
		int a = suf[i-1], b = suf[i], k = 0;
		while (a + k < n && b + k < n && s[a+k] == s[b+k]) k++;
		distinct -= k;
	}
	assert(distinct > INT_MAX);
	assert(sa.distinct() == distinct);
	assert(sa.kth(0) == "a");
	assert(sa.kth(distinct - 1) == s.substr(suf[n-1]));
	assert(sa.kth(distinct - 2) == s.substr(suf[n-1], sz(s) - suf[n-1] - 1));
	assert(sa.kth(1) == (count(all(s), 'a') > 0 && s.find("aa") != string::npos ? "aa" : "ab"));
	vi cnt = sa.occ(s);
	assert(cnt[p] == 1 && cnt[walk(sa, "a")] == count(all(s), 'a'));
	ll tot = 0;
	rep(v,1,sa.N) tot += (ll)cnt[v] * (sa.dis[v] - sa.dis[sa.link[v]]);
	assert(tot == (ll)n * (n + 1) / 2);
	assert(sa.lcs(s) == pii(n, 0));
	assert(sa.lcs("xyz" + s.substr(100, 50) + "xyz") == pii(50, 3));
	string t = randString(n, 3, false);
	pii r = sa.lcs(t);
	assert(walk(sa, t.substr(r.second, r.first)) != -1);
	assert(r.second + r.first == sz(t) || walk(sa, t.substr(r.second, r.first + 1)) == -1);
}

void testPerf() {
	string s = randString(1000000, 26, false);
	SA sa;
	int p = 0;
	for (char c : s) p = sa.append(p, c);
	assert(sa.N <= 2 * sz(s));
	ll edges = 0;
	rep(i,0,sa.N) edges += sz(sa.adj[i]);
	assert(edges <= 3 * sz(s));
	sa.order(); sa.occ(s); sa.distinct(); sa.kth(sa.distinct() - 1); sa.lcs(s);
}

int main() {
	testEdge();
	rep(it,0,6000) {
		int k = 1 + rand() % 3, alpha = 1 + rand() % 3;
		bool bytes = rand() % 10 == 0;
		vector<string> strs, probes;
		rep(i,0,k) {
			int n = rand() % 11;
			string s = randString(n, alpha, bytes);
			if (i && rand() % 3 == 0) {
				string& o = strs[rand() % i];
				int a = rand() % (sz(o) + 1), b = rand() % (sz(o) + 1);
				if (a > b) swap(a, b);
				s = o.substr(a, b - a);
			}
			strs.push_back(s);
		}
		rep(i,0,4) probes.push_back(randString(rand() % 8, alpha + rand() % 2, bytes));
		check(strs, probes);
	}
	testOverflow();
	testPerf();
	cout<<"Tests passed!"<<endl;
}
