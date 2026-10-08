/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: own work; folklore, see https://codeforces.com/blog/entry/45578
 * Description: Operations $0..n-1$ are applied in order to a data
 * structure. For each query $j \in [0,q)$ finds the smallest $k \in [0,n]$
 * such that \texttt{check(j)} holds after the first $k$ operations
 * (i.e. operations $[0,k)$), or $n+1$ if there is none.
 * \texttt{check(j)} must be monotone in $k$ (false, then true) and must
 * not modify the structure. \texttt{reset()} restores the empty structure,
 * \texttt{apply(i)} applies operation $i$.
 * Usage:
	UF uf(0); // edges (a[i], b[i]), queries (u[j], v[j])
	vi r = parBinSearch(m, q, [\&]() { uf = UF(n); },
		[\&](int i) { uf.join(a[i], b[i]); },
		[\&](int j) { return uf.sameSet(u[j], v[j]); });
 * Time: $\lfloor\log_2(n+1)\rfloor+1$ rounds, each one reset, $n$ applies
 * and at most $q$ checks: O((n+q) \log n) operations.
 * Status: stress-tested
 */
#pragma once

template<class R, class A, class C>
vi parBinSearch(int n, int q, R reset, A apply, C check) {
	vi lo(q), hi(q, n + 1);
	vector<vi> b(n + 1);
	for (int r = n + 1; r; r /= 2) {
		rep(j,0,q) if (lo[j] < hi[j])
			b[(lo[j] + hi[j]) / 2].push_back(j);
		reset();
		rep(i,0,n+1) {
			for (int j : b[i])
				if (check(j)) hi[j] = i;
				else lo[j] = i + 1;
			b[i].clear();
			if (i < n) apply(i);
		}
	}
	return lo;
}
