#include "../utilities/template.h"

#include "../../content/strings/AhoCorasickJacob.h"

mt19937 rng(71);

string randStr(int n, int alpha, bool high) {
	string s(n, 'a');
	for (char& c : s) c = (char)((high ? 250 : 'a') + rng() % alpha);
	return s;
}

bool endsWith(const string& s, const string& t) {
	return sz(s) >= sz(t) && s.compare(sz(s) - sz(t), sz(t), t) == 0;
}

void test(const vector<string>& pats, const vector<string>& texts) {
	AhoCorasick ac;
	for (auto& p : pats) ac.add_string(p);
	ac.build();

	// Trie structure: exactly one node per distinct prefix.
	map<string, int> id;
	vector<string> str(sz(ac.nodes));
	function<void(int)> dfs = [&](int v) {
		assert(!id.count(str[v]));
		id[str[v]] = v;
		for (auto p : ac.nodes[v].next) {
			str[p.second] = str[v] + p.first;
			dfs(p.second);
		}
	};
	dfs(0);
	set<string> pref = {""};
	for (auto& p : pats) rep(l,0,sz(p)+1) pref.insert(p.substr(0, l));
	assert(sz(pref) == sz(id) && sz(id) == sz(ac.nodes));
	for (auto& p : pref) assert(id.count(p));

	// Longest suffix of s, of length <= maxLen, that is a node of the trie.
	auto longest = [&](const string& s, int maxLen) {
		for (int l = min(maxLen, sz(s)); l > 0; l--) {
			auto it = id.find(s.substr(sz(s) - l));
			if (it != id.end()) return it->second;
		}
		return 0;
	};
	// link = longest proper suffix that is in the trie
	assert(ac.nodes[0].link == -1);
	rep(v,1,sz(ac.nodes))
		assert(ac.nodes[v].link == longest(str[v], sz(str[v]) - 1));

	set<char> alphabet;
	for (auto& p : pats) for (char c : p) alphabet.insert(c);
	for (auto& t : texts) for (char c : t) alphabet.insert(c);
	rep(v,0,sz(ac.nodes)) for (char c : alphabet)
		assert(ac.nx(v, c) == longest(str[v] + c, sz(str[v]) + 1));

	// match = first index at which some pattern ends
	for (auto& t : texts) {
		int want = -1;
		rep(i,0,sz(t)) {
			string cur = t.substr(0, i + 1);
			for (auto& p : pats) if (endsWith(cur, p)) want = i;
			if (want != -1) break;
		}
		int got = ac.match(t);
		if (got != want) {
			cerr << "match(\"" << t << "\") = " << got << ", expected " << want << "\npatterns:";
			for (auto& p : pats) cerr << " \"" << p << "\"";
			cerr << endl;
			abort();
		}
	}
}

int main() {
	test({}, {"", "a", "abc"});
	test({"a"}, {"", "a", "b", "ba"});
	test({"a", "a"}, {"a", "bbba"});
	test({"abcd"}, {"abc", "abcd", "xabcabcd"});
	// a pattern that is a proper suffix of a prefix of another one
	test({"abc", "b"}, {"ab", "abx", "abc", "xxb"});
	test({"aaaa", "ab", "bab", "b"}, {"aaab", "aaaa", "cccc"});

	rep(it,0,60000) {
		int alpha = (int)(rng() % 3) + 1;
		bool high = rng() % 8 == 0;
		int k = (int)(rng() % 5);
		vector<string> pats, texts;
		rep(i,0,k) pats.push_back(randStr((int)(rng() % 5) + 1, alpha, high));
		if (k && rng() % 4 == 0) pats.push_back(pats[0]); // duplicate
		rep(i,0,4) texts.push_back(randStr((int)(rng() % 12), alpha + (int)(rng() % 2), high));
		test(pats, texts);
	}
	// longer patterns, so that most texts do not match early
	rep(it,0,3000) {
		int k = (int)(rng() % 6) + 1;
		vector<string> pats, texts;
		rep(i,0,k) pats.push_back(randStr((int)(rng() % 8) + 3, 2, false));
		rep(i,0,4) texts.push_back(randStr((int)(rng() % 40), 2, false));
		test(pats, texts);
	}
	cout<<"Tests passed!"<<endl;
}
