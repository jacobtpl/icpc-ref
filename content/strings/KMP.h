/**
 * Author: Johan Sannemo
 * Date: 2016-12-15
 * License: CC0
 * Description: pi[x] computes the length of the longest prefix of s that ends at x, other than s[0...x] itself (abacaba -> 0010123).
 * Can be used to find all occurrences of a string.
 * Time: O(n)
 * Status: Tested on kattis:stringmatching
 */
#pragma once

vi pi(const string& s) {
	vi p(sz(s));
	rep(i,1,sz(s)) {
		int g = p[i-1];
		while (g && s[i] != s[g]) g = p[g-1];
		p[i] = g + (s[i] == s[g]);
	}
	return p;
}
vi match(const string& s, const string& pat) {
	vi p = pi(pat + '\0' + s), res;
	rep(i,sz(p)-sz(s),sz(p))
		if (p[i] == sz(pat)) res.push_back(i - 2 * sz(pat));
	return res;
}
vi match2(const string& s, const string& pat) { // only compute pi for pat
	vi p = pi(pat), res;
	int cp = 0;
	rep(i,1,sz(s)) {
		int g = cp;
		while (g && s[i] != pat[g]) g = p[g-1];
		cp = g + (s[i] == pat[g]);
		if (cp >= sz(pat)) res.push_back(i - sz(pat) + 1);
	}
	return res;
}