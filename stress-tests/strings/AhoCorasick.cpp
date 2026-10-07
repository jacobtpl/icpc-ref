#include "../utilities/template.h"
#include "../../content/strings/AhoCorasick.h"

#define trav(a, x) for (auto& a : x)

template<class F>
void gen(string& s, int at, int alpha, F f) {
	if (at == sz(s)) f();
	else {
		rep(i,0,alpha) {
			s[at] = (char)('A' + i);
			gen(s, at+1, alpha, f);
		}
	}
}

void test(const string& s) {
	vector<string> pats;
	string cur;
	rep(i,0,sz(s)) {
		if (s[i] == 'A') {
			pats.push_back(cur);
			cur = "";
		}
		else cur += s[i];
	}

	string hay = cur;
	trav(x, pats) if (x.empty()) return;

	AhoCorasick ac(pats);
	vector<vi> positions = ac.findAll(pats, hay);

	vi ord;
	rep(i,0,sz(hay)) {
		ord.clear();
		rep(j,0,sz(pats)) {
			string& pat = pats[j];
			if (hay.substr(i, pat.size()) == pat) {
				ord.push_back(j);
			}
		}
		sort(all(positions[i]));

		if (positions[i] != ord) {
			cerr << "failed!" << endl;
			cerr << hay << endl;
			trav(x, pats) cerr << x << endl;
			cerr << "failed at position " << i << endl;
			cerr << "got:" << endl;
			trav(x, positions[i]) cerr << x << ' ';
			cerr << endl;
			cerr << "expected:" << endl;
			trav(x, ord) cerr << x << ' ';
			cerr << endl;
			abort();
		}
	}
}

// find() and nmatches against the definition, on random inputs
// (duplicate patterns, patterns that are suffixes/prefixes of each other).
void testRandom() {
	mt19937 rng(4242);
	rep(it,0,100000) {
		int alpha = (int)(rng() % 3) + 1, k = (int)(rng() % 6);
		if (it % 100 == 0) alpha = 26;
		vector<string> pats(k);
		for (auto& p : pats) {
			p.resize(rng() % 5 + 1);
			for (char& c : p) c = (char)('A' + (alpha == 26 && rng() % 2 ? 25 : rng() % alpha));
		}
		if (k && rng() % 3 == 0) pats.push_back(pats[rng() % k]);
		string hay(rng() % 20, 'A');
		for (char& c : hay) c = (char)('A' + (alpha == 26 && rng() % 2 ? 25 : rng() % alpha));

		AhoCorasick ac(pats);
		vi r = ac.find(hay);
		vector<vi> all = ac.findAll(pats, hay);
		assert(sz(r) == sz(hay) && sz(all) == sz(hay));
		int n = 0;
		rep(i,0,sz(hay)) {
			int longest = -1, cnt = 0;
			rep(j,0,sz(pats)) {
				int l = sz(pats[j]);
				if (l <= i + 1 && hay.compare(i + 1 - l, l, pats[j]) == 0) {
					cnt++;
					if (longest == -1 || l > sz(pats[longest])) longest = j;
				}
			}
			if (longest == -1) assert(r[i] == -1);
			else assert(r[i] != -1 && pats[r[i]] == pats[longest]);
			n = ac.N[n].next[hay[i] - AhoCorasick::first];
			assert(ac.N[n].nmatches == cnt);

			vi ord;
			rep(j,0,sz(pats)) if (hay.compare(i, sz(pats[j]), pats[j]) == 0) ord.push_back(j);
			// "shortest first"
			rep(j,1,sz(all[i])) assert(sz(pats[all[i][j-1]]) <= sz(pats[all[i][j]]));
			sort(all(all[i]));
			assert(all[i] == ord);
		}
	}
}

int main() {
	testRandom();
	// test ~4^10 strings
	rep(n,0,11) {
		string s(n, 'x');
		gen(s, 0, 4, [&]() {
			test(s);
		});
	}
	// then ~5^7
	rep(n,0,8) {
		string s(n, 'x');
		gen(s, 0, 5, [&]() {
			test(s);
		});
	}
	cout<<"Tests passed!"<<endl;
}
