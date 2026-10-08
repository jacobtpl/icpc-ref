/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: Duval (1983), https://cp-algorithms.com/string/lyndon_factorization.html
 * Description: Lyndon factorisation (Duval): the unique split
 *  $s = w_1 w_2 \dots w_k$ with $w_1 \ge w_2 \ge \dots \ge w_k$, each $w_i$
 *  strictly smaller than all its proper suffixes.
 *  Returns the $k+1$ boundaries: 0-indexed, $w_i$ = s[r[i-1], r[i]),
 *  r[0] = 0, r[k] = n (just \{0\} for empty $s$).
 *  S is a string or vector, compared with operator<.
 *  r[k-1] is the start of the smallest nonempty suffix of $s$.
 *  Min rotation of $s$ starts at the last r[i] $< n$ in lyndon(s+s).
 * Time: O(n)
 * Usage:
 *  vi r = lyndon(s); // "abaab" -> 0 2 5 (ab|aab)
 * Status: stress-tested
 */
#pragma once

template<class S> vi lyndon(const S& s) {
	int n = sz(s), i = 0;
	vi r;
	while (i < n) {
		int j = i + 1, k = i;
		for (; j < n && !(s[j] < s[k]); j++)
			k = s[k] < s[j] ? i : k + 1;
		for (; i <= k; i += j - k) r.push_back(i);
	}
	r.push_back(n);
	return r;
}
