#include "../utilities/template.h"

#include "../../content/strings/SuffixAutomaton.h"

set<string> substrings(const vector<string>& v) {
	set<string> r;
	for (auto& s : v) rep(i,0,sz(s)) rep(l,1,sz(s)-i+1)
		r.insert(s.substr(i, l));
	return r;
}

// every path from the root, as (string, node)
void paths(SA& sa, int v, string& cur, vector<pair<string,int>>& out) {
	for (auto [c, u] : sa.adj[v]) {
		cur += c;
		out.push_back({cur, u});
		paths(sa, u, cur, out);
		cur.pop_back();
	}
}

void test(const vector<string>& v) {
	SA sa;
	int total = 0, nonEmpty = 0;
	for (auto& s : v) {
		int p = 0;
		for (char c : s) {
			p = sa.append(p, c);
			assert(0 <= p && p < sa.N);
		}
		assert(sa.dis[p] == sz(s));
		total += sz(s), nonEmpty += !s.empty();
	}
	assert(sa.N == sz(sa.adj) && sa.N == sz(sa.link) && sa.N == sz(sa.dis));
	assert(sa.N <= max(1, 2 * total));
	assert(sa.link[0] == -1 && sa.dis[0] == 0);
	set<string> subs = substrings(v);
	// distinct substrings = sum of dis[v] - dis[link[v]]
	ll cnt = 0;
	rep(i,1,sa.N) {
		assert(0 <= sa.link[i] && sa.link[i] < sa.N);
		assert(sa.dis[sa.link[i]] < sa.dis[i]);
		cnt += sa.dis[i] - sa.dis[sa.link[i]];
	}
	assert(cnt == sz(subs));
	// the language of the automaton is exactly the set of substrings
	vector<pair<string,int>> all;
	string cur;
	paths(sa, 0, cur, all);
	assert(sz(all) == sz(subs));
	vi mx(sa.N, 0), mn(sa.N, INT_MAX), num(sa.N);
	vector<string> longest(sa.N);
	for (auto& [s, u] : all) {
		assert(subs.count(s));
		num[u]++;
		mn[u] = min(mn[u], sz(s));
		if (sz(s) > mx[u]) mx[u] = sz(s), longest[u] = s;
	}
	rep(i,1,sa.N) if (num[i]) {
		// a state holds the suffixes of its longest string with
		// lengths in (dis[link], dis]
		assert(mx[i] == sa.dis[i]);
		assert(mn[i] == sa.dis[sa.link[i]] + 1);
		assert(num[i] == mx[i] - mn[i] + 1);
		int l = sa.link[i];
		if (l) assert(longest[l] ==
			longest[i].substr(sz(longest[i]) - sa.dis[l]));
	}
	for (auto& [s, u] : all)
		assert(longest[u].compare(sz(longest[u]) - sz(s), sz(s), s) == 0);
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
	test({});
	test({""});
	test({"", "", "a", ""});
	test({"a", "a", "a"});
	// exhaustive single strings
	rep(n,0,13) { string s(n, 'x'); gen(s, 0, 2, [&]() { test({s}); }); }
	rep(n,0,9) { string s(n, 'x'); gen(s, 0, 3, [&]() { test({s}); }); }
	// exhaustive pairs (generalized automaton)
	rep(n,0,7) rep(m,0,7) {
		string s(n, 'x'), t(m, 'x');
		gen(s, 0, 2, [&]() { gen(t, 0, 2, [&]() { test({s, t}); }); });
	}
	mt19937 rng(4242);
	auto rnd = [&](int a, int b) {
		return (int)(rng() % (unsigned)(b - a + 1)) + a; };
	rep(it,0,30000) {
		int k = rnd(1, 4), alpha = rnd(1, 4), big = it % 100 == 0;
		vector<string> v(k);
		for (auto& s : v) {
			int n = rnd(0, big ? 60 : 12);
			rep(i,0,n) s += (char)('a' + rnd(0, alpha-1));
			if (rnd(0, 4) == 0) {
				int p = rnd(1, 4);
				rep(i,p,n) s[i] = s[i-p];
			}
		}
		if (k > 1 && rnd(0, 3) == 0) v[1] = v[0]; // duplicates
		if (k > 2 && rnd(0, 3) == 0) // prefix / suffix of another
			v[2] = rnd(0, 1) ? v[0].substr(0, sz(v[0]) / 2)
				: v[0].substr(sz(v[0]) / 2);
		test(v);
	}
	// non-letter bytes are fine as keys of map<char,int>
	test({string("\x01\xff\x80\x01\xff", 5), string("\xff\x80\x7f", 3)});
	cout<<"Tests passed!"<<endl;
}
