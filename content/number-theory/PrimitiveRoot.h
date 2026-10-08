/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: own; https://cp-algorithms.com/algebra/primitive-root.html
 * Description: Smallest primitive root $g$ modulo a prime $p < 7 \cdot 10^{18}$:
 * $g$ is one iff $g^{(p-1)/q} \neq 1$ for every prime $q \mid p-1$.
 * For odd $p$, $h = g$ (or $g+p$ if $g^{p-1} \equiv 1 \pmod{p^2}$) is a
 * primitive root mod every $p^k$, and the odd one of $h$, $h+p^k$ mod $2p^k$.
 * Time: $O(p^{1/4})$ to factor $p-1$, then $O(\log^2 p)$ per candidate
 * ($g \le 113$ for all $p < 10^8$).
 * Status: stress-tested
 */
#pragma once

#include "Factor.h"

ull primRoot(ull p) {
	if (p == 2) return 1;
	auto f = factor(p - 1);
	for (ull g = 2;; g++) {
		bool ok = 1;
		for (ull q : f) if (modpow(g, (p-1) / q, p) == 1) ok = 0;
		if (ok) return g;
	}
}
