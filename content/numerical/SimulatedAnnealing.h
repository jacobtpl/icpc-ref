/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: Kirkpatrick, Gelatt, Vecchi (1983); splitmix64 (Steele, Lea, Flood 2014; Vigna's constants)
 * Description: Simulated annealing, minimizing a score for \texttt{sec} seconds.
 * \texttt{step()} applies a random neighbour move to the user's state and returns
 * the new score, \texttt{undo()} reverts that move, \texttt{save()} is called
 * whenever the state is strictly better than all earlier ones (copy it there); the
 * initial state is not saved. Returns the best score, given the initial score \texttt{cur}.
 * A worse move is accepted with probability $e^{-\Delta/T}$, where $T$ falls
 * geometrically from \texttt{t0} to \texttt{t1} ($\texttt{t0} \ge \texttt{t1} > 0$);
 * pick them near the largest and smallest typical $\Delta$ of a worsening move.
 * Scores must be finite. To maximize, negate the score.
 * The clock is read every 64 steps, so it overruns by up to 64 steps: keep step fast.
 * \texttt{rnd()} is uniform on $[0, 2^{64})$ with all bits usable (\texttt{rnd() \% n} is fine
 * for any $n$), \texttt{rnd01()} on $[0, 1)$; any seed \texttt{rs} works.
 * Time: \texttt{sec} seconds
 * Usage:
	double b = anneal(len(p), 0.9, 100, 0.1, [\&]() {
		j = rnd() \% n, i = rnd() \% (j + 1); // 2-opt, say
		rev(i, j); return len(p); // better: update in O(1)
	}, [\&]() { rev(i, j); }, [\&]() { bp = p; });
 * Status: stress-tested (optimal on small TSP against DP, Boltzmann occupancy, RNG pair coverage and chi-square)
 */
#pragma once

uint64_t rs = 88172645463325252;
uint64_t rnd() {
	auto z = rs += 0x9e3779b97f4a7c15;
	z = (z ^ z >> 30) * 0xbf58476d1ce4e5b9;
	z = (z ^ z >> 27) * 0x94d049bb133111eb;
	return z ^ z >> 31;
}
double rnd01() { return double(rnd() >> 11) / (1LL << 53); }

template<class S, class U, class B>
double anneal(double cur, double sec, double t0, double t1,
		S step, U undo, B save) {
	auto st = chrono::steady_clock::now();
	double best = cur, T = t0;
	for (ll it = 0;; it++) {
		if (!(it & 63)) {
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
