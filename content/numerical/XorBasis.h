/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: own work; folklore (linear basis over $\mathbb F_2$)
 * Description: Basis of a set of 64-bit vectors over $\mathbb F_2$.
 *  $b[i]$ is 0 or a basis vector with highest set bit $i$; $r$ is the rank,
 *  so $2^r$ distinct values are representable as xors of subsets (0 always is).
 *  add returns true iff $x$ was not representable before.
 *  maxXor/minXor give the max/min of $x \oplus v$ over representable $v$.
 *  kth returns the $k$-th smallest representable value, 0-indexed
 *  ($\texttt{kth(0)} = 0$), and requires $k < 2^r$.
 *  Values are compared as unsigned; cast negative numbers to ull.
 * Time: O(64) per operation, O(64^2) for merge.
 * Usage:
 *  XorBasis B; B.add(5); B.add(3);
 *  B.has(6); // true
 *  B.maxXor(); B.kth(2); // 6, 5
 * Status: stress-tested against subset enumeration
 */
#pragma once

typedef unsigned long long ull;
struct XorBasis {
	ull b[64] = {}; int r = 0;
	bool add(ull x) {
		for (int i = 64; i--;) if (x >> i & 1) {
			if (!b[i]) return b[i] = x, ++r;
			x ^= b[i];
		}
		return 0;
	}
	ull minXor(ull x) {
		for (int i = 64; i--;) x = min(x, x ^ b[i]);
		return x;
	}
	ull maxXor(ull x = 0) {
		for (int i = 64; i--;) x = max(x, x ^ b[i]);
		return x;
	}
	bool has(ull x) { return !minXor(x); }
	ull kth(ull k) {
		ull x = 0; int c = r;
		for (int i = 64; i--;) if (b[i])
			if ((x >> i ^ k >> --c) & 1) x ^= b[i];
		return x;
	}
	void merge(const XorBasis& o) { for (ull x : o.b) add(x); }
};
