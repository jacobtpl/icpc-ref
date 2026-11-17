/**
 * Author: chilli, Ramchandra Apte, Noam527, Simon Lindholm
 * Date: 2019-04-24
 * License: CC0
 * Source: https://github.com/RamchandraApte/OmniTemplate/blob/master/modulo.hpp
 * Description: Calculate a* bmod c (or a^b mod c) for 0 <= a, b <= c <= 7.2* 10^18.
 * Time: O(1) for modmul, O(log b) for modpow
 * Status: stress-tested, proven correct
 * Details:
 * This runs ~2x faster than the naive (__int128_t)a * b % M.
 * A proof of correctness is in doc/modmul-proof.tex. An earlier version of the proof,
 * from when the code used a * b / (long double)M, is in doc/modmul-proof.md.
 * The proof assumes that long doubles are x87 80-bit floats (x86-64 Linux). Where they
 * are not (e.g. Apple Silicon, MSVC), the __int128 version is used instead.
 */
#pragma once

typedef unsigned long long ull;
ull modmul(ull a, ull b, ull M) {
#if __LDBL_MANT_DIG__ != 64
	return ull((__uint128_t)a * b % M);
#endif
	ll ret = a * b - M * ull(1.L / M * a * b);
	return ret + M * (ret < 0) - M * (ret >= (ll)M);
}
ull modpow(ull b, ull e, ull mod) {
	ull ans = 1;
	for (; e; b = modmul(b, b, mod), e /= 2)
		if (e & 1) ans = modmul(ans, b, mod);
	return ans;
}
