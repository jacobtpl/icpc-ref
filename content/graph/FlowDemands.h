/**
 * Author: Devin
 * Date: 2026-10-08
 * License: CC0
 * Source: own work; standard reduction (super source/sink for the
 *  lower bounds, then augmenting $s \to t$ in the residual graph)
 * Description: Flow with edge bounds $lo \le f \le hi$ ($lo$ may be negative).
 *  Nodes are 0-indexed; \texttt{addEdge} returns the edge's index $i$
 *  and \texttt{flow(i)} is its flow after a successful solve.
 *  \texttt{circulation()} tells whether a feasible circulation exists.
 *  \texttt{maxFlow(s,t)}, $s \ne t$, returns the maximum value (net flow out of $s$)
 *  of a feasible $s$-$t$ flow of value $\ge 0$, or $-1$ if there is none;
 *  it must be bounded. Node supplies: do \texttt{ex[v] += d} before solving
 *  to make $out(v) - in(v) = d$ (they must sum to $0$).
 *  Call only one solve, once. Sums of $|lo|$ and of $hi-lo$ must fit in ll.
 * Time: one or two runs of Dinic ($O(V^2E)$) on $V+2$ nodes and $E+V+1$ edges.
 * Usage: FlowDemands F(n); int i = F.addEdge(u, v, lo, hi);
 *  if (F.maxFlow(s, t) >= 0) x = F.flow(i);
 * Status: stress-tested against brute force and cut conditions
 */
#pragma once

#include "Dinic.h"

struct FlowDemands {
	int n; Dinic D;
	vector<ll> ex, lo; vi id;
	FlowDemands(int n) : n(n), D(n + 2), ex(n) {}
	int addEdge(int u, int v, ll l, ll h) {
		id.pb(u == v ? -1 : sz(D.E)); lo.pb(l);
		D.AddEdge(u, v, h - l);
		ex[v] += l; ex[u] -= l;
		return sz(lo) - 1;
	}
	bool circulation() {
		ll need = 0;
		rep(i,0,n)
			if (ex[i] > 0) D.AddEdge(n, i, ex[i]), need += ex[i];
			else D.AddEdge(i, n + 1, -ex[i]);
		return D.MaxFlow(n, n + 1) == need;
	}
	ll maxFlow(int s, int t) {
		addEdge(t, s, 0, LLONG_MAX);
		return circulation() ? D.MaxFlow(s, t) : -1;
	}
	ll flow(int i) { return lo[i] + (id[i] < 0 ? 0 : D.E[id[i]].flow); }
};
