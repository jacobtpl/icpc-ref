/**
 * Author: Devin
 * Date: 2026-11-22
 * License: CC0
 * Source: UCup Season 5 Stage 8 E; Idury-Schaffer parameterized matching
 * Description: Aho-Corasick up to renaming of letters. A pattern p matches
 * a window w if there is a bijection f on letters with f(p[i]) = w[i]
 * (same equality pattern). Works on any string/vector with comparable items.
 * add(p, w) adds weight w to pattern p; call build() once after all adds.
 * After build, t[v].val = total weight of patterns matching a suffix of the
 * text read so far. run(s, f) calls f(i, v) with the node after reading s[i].
 * For some letters fixed (must match exactly), give them distance 0 in dist
 * and add their identity to the code.
 * Time: O((sum |p| + |s|) * K), K = distinct letters + 1 (child list scan).
 * Status: stress-tested
 */
#pragma once

struct ParamAC {
	struct Node { int len, c, sib, head = -1, link = 0; ll val = 0; };
	vector<Node> t{{0, 0, -1}};
	template<class S> static vi dist(const S& s) { // i - prev equal, or 0
		map<typename S::value_type, int> last; vi d(sz(s));
		rep(i,0,sz(s)) {
			auto it = last.find(s[i]);
			if (it != last.end()) d[i] = i - it->second, it->second = i;
			else last[s[i]] = i;
		}
		return d;
	}
	int child(int v, int c) {
		for (int x = t[v].head; x != -1; x = t[x].sib)
			if (t[x].c == c) return x;
		return -1;
	}
	template<class S> int add(const S& s, ll w = 1) {
		vi d = dist(s); int v = 0;
		rep(i,0,sz(s)) {
			int x = child(v, d[i]);
			if (x < 0) {
				x = sz(t); t.push_back({i + 1, d[i], t[v].head});
				t[v].head = x;
			}
			v = x;
		}
		t[v].val += w;
		return v;
	}
	int step(int v, int d) { // d = distance to previous equal item, 0 if none
		for (;; v = t[v].link) {
			int x = child(v, d > t[v].len ? 0 : d);
			if (x >= 0) return x;
			if (!v) return 0;
		}
	}
	void build() {
		vi q{0};
		rep(i,0,sz(q)) {
			int v = q[i];
			for (int x = t[v].head; x != -1; x = t[x].sib) {
				t[x].link = v ? step(t[v].link, t[x].c) : 0;
				t[x].val += t[t[x].link].val; // combine along suffix links
				q.push_back(x);
			}
		}
	}
	template<class S, class F> void run(const S& s, F f) {
		vi d = dist(s); int v = 0;
		rep(i,0,sz(s)) f(i, v = step(v, d[i]));
	}
};
