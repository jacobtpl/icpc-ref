/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: Kirkpatrick, Gelatt, Vecchi (1983); xorshift64 (7, 9) by Marsaglia
 * Description: Simulated annealing, minimizing a score for \texttt{sec} seconds.
 * \texttt{step()} applies a random neighbour move to the user's state and returns
 * the new score, \texttt{undo()} reverts that move, \texttt{save()} is called
 * whenever the state is the best seen so far (copy it there); the initial state is
 * not saved. Returns the best score, given the initial score \texttt{cur}.
 * A worse move is accepted with probability $e^{-\Delta/T}$, where $T$ falls
 * geometrically from \texttt{t0} to \texttt{t1} ($\texttt{t0} \ge \texttt{t1} > 0$);
 * pick them near the largest and smallest typical $\Delta$ of a worsening move.
 * Scores must be finite. To maximize, negate the score.
 * The clock is read every 256 steps, so step should be fast ($O(1)$ score updates).
 * \texttt{rnd()} is uniform on $[1, 2^{64})$, \texttt{rnd01()} on $[0, 1)$; \texttt{rs} must be nonzero.
 * Time: \texttt{sec} seconds; one step, and one undo or save, per iteration.
 * Usage:
	double b = anneal(len(p), 0.9, 100, 0.1, [&]() {
		i = rnd() % n, j = rnd() % n; // 2-opt, say
		rev(i, j); return len(p); // better: update in O(1)
	}, [&]() { rev(i, j); }, [&]() { bp = p; });
 * Status: stress-tested (optimal on small TSP against DP, Boltzmann acceptance, RNG period)
 */
#pragma once

typedef unsigned long long ull;
ull rs = 88172645463325252;
ull rnd() { rs ^= rs << 7; return rs ^= rs >> 9; }
double rnd01() { return double(rnd() >> 11) * 0x1p-53; }

template<class S, class U, class B>
double anneal(double cur, double sec, double t0, double t1,
		S step, U undo, B save) {
	auto st = chrono::steady_clock::now();
	double best = cur, T = t0;
	for (ull it = 0;; it++) {
		if (!(it & 255)) {
			double t = chrono::duration<double>(
				chrono::steady_clock::now() - st).count();
			if (t >= sec) break;
			T = t0 * pow(t1 / t0, t / sec);
		}
		double nw = step();
		if (nw <= cur || rnd01() < exp((cur - nw) / T)) {
			cur = nw;
			if (cur < best) best = cur, save();
		} else undo();
	}
	return best;
}
