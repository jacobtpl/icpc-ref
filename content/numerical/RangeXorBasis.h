/**
 * Author: Devin
 * Date: 2026-11-17
 * License: CC0
 * Source: folklore (prefix linear basis, CF 1100F); UCup 5 Stage 1 F
 * Description: Xor basis of any subarray a[l..r]. For every prefix it keeps
 *  the basis that prefers the latest positions: b[r][k] is a vector with
 *  highest bit k (or 0) and p[r][k] the position it came from. The vectors
 *  with p[r][k] >= l form a basis of a[l..r] (0-indexed, inclusive).
 *  Values must be < 2^B.
 *  Offline (queries sorted by r) keep only the current row to save memory.
 *  Same idea works for vectors mod a prime: on a pivot clash keep the later
 *  vector, continue eliminating with the earlier one.
 * Memory: 12 * B bytes per element
 * Time: O(B) per push and per query
 * Usage: RangeXorBasis rb; for (ll x : a) rb.push(x);
 *  rb.maxXor(l, r); rb.rank(l, r); rb.has(l, r, x);
 * Status: stress-tested against brute force
 */
#pragma once

struct RangeXorBasis {
	static const int B = 30;
	typedef unsigned T; // ull for B > 32
	vector<array<T, B>> b; vector<array<int, B>> p;
	void push(T x) {
		int i = sz(b);
		b.push_back(i ? b[i-1] : array<T, B>{});
		p.push_back(i ? p[i-1] : array<int, B>{});
		auto &nb = b[i]; auto &np = p[i];
		for (int k = B; k-- && x;) if (x >> k & 1) {
			if (!nb[k]) { nb[k] = x; np[k] = i; break; }
			if (np[k] < i) swap(nb[k], x), swap(np[k], i);
			x ^= nb[k];
		}
	}
	T minXor(int l, int r, T x) { // min of x ^ v, v in span(a[l..r])
		for (int k = B; k--;)
			if (b[r][k] && p[r][k] >= l) x = min(x, x ^ b[r][k]);
		return x;
	}
	T maxXor(int l, int r, T x = 0) {
		for (int k = B; k--;)
			if (b[r][k] && p[r][k] >= l) x = max(x, x ^ b[r][k]);
		return x;
	}
	bool has(int l, int r, T x) { return !minXor(l, r, x); }
	int rank(int l, int r) {
		int c = 0;
		rep(k,0,B) c += b[r][k] && p[r][k] >= l;
		return c;
	}
};
