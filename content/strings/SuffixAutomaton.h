/**
 * Author: Siyong Huang
 * Date: 2023-05-01
 * License: CC0
 * Source: My head
 * Description: Generalized suffix automaton. State 0 is the
 * empty string; \texttt{append(p, c)} returns the state of
 * (a string of state $p$) + $c$. Each state $v \ne 0$ holds the
 * substrings that are suffixes of its longest one with length in
 * $(dis[link[v]], dis[v]]$; they all have the same end
 * positions. At most $2n$ states and $3n$ transitions.
 * See SuffixAutomatonExamples.h for common uses.
 * Time: \texttt{append} is amortized $O(\log \Sigma)$, so
 * building is $O(n \log \Sigma)$.
 * Usage:
 *  SA sa; int p = 0;
 *  for (char c : s) p = sa.append(p, c);
 *  // for another string, start again from p = 0
 * Status: stress-tested
 */
#pragma once

struct SA {
	vector<map<char, int> > adj;
	vi link, dis;
	int N;
	SA(): adj(1), link(1, -1), dis(1, 0), N(1) {}
	int new_node(int v=-1) {
		if(v == -1)
			adj.emplace_back(), link.emplace_back(), dis.emplace_back();
		else
			adj.push_back(adj[v]), link.push_back(link[v]), dis.push_back(dis[v]);
		return N++;
	}
	int go(int p, char c) {
		auto it = adj[p].find(c);
		if(dis[it->second] == dis[p] + 1)
			return it->second;
		else {
			int q = it->second, n = new_node(q);
			dis[n] = dis[p] + 1, link[q] = n;
			for(;p != -1 && (it = adj[p].find(c))->second == q;p = link[p])
				it->second = n;
			return n;
		}
	}
	int append(int p, char c) {
		auto it = adj[p].find(c);
		if(it != adj[p].end())
			return go(p, c);
		int n = new_node();
		dis[n] = dis[p] + 1;
		for(;p != -1 && adj[p].find(c) == adj[p].end();p = link[p])
			adj[p].insert({c, n});
		if(p == -1)
			link[n] = 0;
		else
			link[n] = go(p, c);
		return n;
	}
};
