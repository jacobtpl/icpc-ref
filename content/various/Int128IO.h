/**
 * Author: Devin
 * License: CC0
 * Description: Decimal conversion for the full signed 128-bit range.
 *  parse128 accepts an optional sign followed by digits. It returns false
 *  on invalid input or overflow, leaving out unchanged. str128 handles
 *  the minimum value without signed negation overflow.
 * Time: O(|s|) parsing, O(\log |x|) formatting.
 * Status: Stress-tested against Python integers.
 */
#pragma once

bool parse128(const string& s, __int128_t& out) {
	if (s.empty()) return false;
	bool neg = s[0] == '-';
	size_t i = neg || s[0] == '+';
	if (i == s.size()) return false;
	using U = __uint128_t;
	U n = 0, limit = (U(1) << 127) - !neg;
	for (; i < s.size(); ++i) {
		unsigned d = unsigned(s[i] - '0');
		if (d > 9 || n > (limit-d)/10) return false;
		n = n*10 + d;
	}
	out = neg && n ? -__int128_t(n-1)-1 : __int128_t(n);
	return true;
}

string str128(__int128_t x) {
	using U = __uint128_t;
	U n = x < 0 ? -U(x) : U(x);
	string s;
	do {
		s += char('0' + n%10);
		n /= 10;
	} while (n);
	if (x < 0) s += '-';
	reverse(all(s));
	return s;
}
