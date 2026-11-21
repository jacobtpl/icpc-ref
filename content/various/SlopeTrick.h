/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: own work; idea from zscoder's Codeforces blog "Slope trick
 *  explained"
 * Description: Maintains a convex piecewise linear function $f$
 *  (over integers or reals), initially $f = 0$, storing the
 *  breakpoints left/right of the minimum in two heaps.
 *  \texttt{mn} is $\min f$, attained exactly on
 *  $[\texttt{topL()}, \texttt{topR()}]$ ($\mp$\texttt{inf} if unbounded).
 *  \texttt{shift(a, b)} requires $a \le b$; with $a = b$ it
 *  translates $f$ right by $a$. The largest $|a|$ plus the sum of
 *  all $|$shift arguments$|$ must be $\le 10^{18}$, and \texttt{mn}
 *  must fit in a \texttt{ll}.
 * Time: $O(\log N)$ per add, $O(1)$ for shift,
 *  $O(N)$ for prefMin/sufMin (amortized $O(1)$).
 * Usage:
 *  SlopeTrick f; // min sum |a_i-b_i|, b non-decreasing:
 *  for (ll x : a) f.prefMin(), f.addAbs(x);
 *  ll cost = f.mn;
 * Status: stress-tested
 */
#pragma once

struct SlopeTrick {
	static constexpr ll inf = LLONG_MAX / 4;
	ll mn = 0, aL = 0, aR = 0;
	priority_queue<ll> L;
	priority_queue<ll, vector<ll>, greater<ll>> R;
	ll topL() { return sz(L) ? L.top() + aL : -inf; }
	ll topR() { return sz(R) ? R.top() + aR : inf; }
	void add(ll c) { mn += c; } // f += c
	void addR(ll a) { // f += max(x - a, 0)
		mn += max(0LL, topL() - a);
		L.push(a - aL); R.push(topL() - aR); L.pop();
	}
	void addL(ll a) { // f += max(a - x, 0)
		mn += max(0LL, a - topR());
		R.push(a - aR); L.push(topR() - aL); R.pop();
	}
	void addAbs(ll a) { addR(a); addL(a); } // f += |x - a|
	void prefMin() { R = {}; } // f(x) = min f(y), y <= x
	void sufMin() { L = {}; } // f(x) = min f(y), y >= x
	// f(x) = min f(y), x - b <= y <= x - a
	void shift(ll a, ll b) { aL += a; aR += b; }
};
