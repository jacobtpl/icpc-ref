/**
 * Author: Siyong Huang
 * Date: 2023-05-01
 * License: CC0
 * Source: My head
 * Description: Generalized suffix automaton (several strings may
 * be added) over characters of type T (SA<int> for int alphabets).
 * State 0 is the empty string. Each state v != 0 holds
 * the substrings that are suffixes of its longest one with length
 * in (dis[link[v]], dis[v]]; they all have the same end positions.
 * link[v] is the state of the longest suffix not in v (link[0] = -1).
 * adj[v][c] is the state of (string of v) + c. At most 2n states
 * and 3n transitions. See SuffixAutomatonExamples.h for how to use it.
 * Time: append is amortized O(log sigma), building is O(n log sigma).
 * Usage:
 *  SA sa; int p = 0; // SA<int> sa; for basic_string<int>
 *  for (char c : s) p = sa.append(p, c); // p = state of s
 *  p = 0; for (char c : t) p = sa.append(p, c); // add another string
 *  // is w a substring of some added string?
 *  int v = 0;
 *  for (char c : w) { auto it = sa.adj[v].find(c);
 *    if (it == sa.adj[v].end()) { v = -1; break; } v = it->second; }
 *  bool found = v != -1;
 * Status: stress-tested
 */
#pragma once

template<class T = char> struct SA {
	vector<map<T, int>> adj;
	vi link, dis;
	int N = 1;
	SA(): adj(1), link(1, -1), dis(1, 0) {}
	int node(int v = -1) {
		if (v == -1) adj.emplace_back(), link.push_back(0), dis.push_back(0);
		else adj.push_back(adj[v]), link.push_back(link[v]), dis.push_back(dis[v]);
		return N++;
	}
	int go(int p, T c, typename map<T, int>::iterator it) {
		int q = it->second;
		if (dis[q] == dis[p] + 1) return q;
		int n = node(q);
		dis[n] = dis[p] + 1, link[q] = n;
		do it->second = n;
		while ((p = link[p]) != -1 && (it = adj[p].find(c))->second == q);
		return n;
	}
	int append(int p, T c) {
		auto it = adj[p].find(c);
		if (it != adj[p].end()) return go(p, c, it);
		int n = node();
		dis[n] = dis[p] + 1;
		while (p != -1 && (it = adj[p].emplace(c, n).first)->second == n)
			p = link[p];
		link[n] = p == -1 ? 0 : go(p, c, it);
		return n;
	}
};
