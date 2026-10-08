/**
 * Author: Siyong Huang, Devin
 * Date: 2023-05-01
 * License: CC0
 * Source: My head; extras after
 *  https://cp-algorithms.com/string/suffix-automaton.html
 * Description: Generalized suffix automaton. State 0 is the
 * empty string; \texttt{append(p, c)} returns the state of
 * (a string of state $p$) + $c$. Each state $v \ne 0$ holds the
 * substrings that are suffixes of its longest one with length in
 * $(dis[link[v]], dis[v]]$; they all have the same end
 * positions. At most $2n$ states and $3n$ transitions.
 * \texttt{order}: states by increasing \texttt{dis} (a
 * topological order of \texttt{adj}, children before parents
 * in the link tree when reversed).
 * \texttt{occ(s)}: for every state, the number of occurrences
 * in $s$ of each of its substrings ($s$ must have been added;
 * entry 0 is meaningless).
 * \texttt{distinct}: number of distinct non-empty substrings.
 * \texttt{kth(k)}: $k$-th (0-indexed, $0 \le k <$
 * \texttt{distinct()}) distinct non-empty substring in sorted
 * order (by \texttt{char}, which is usually signed).
 * \texttt{lcs(t)}: \{length, start in $t$\} of a longest
 * substring of $t$ that is a substring of an added string.
 * Time: \texttt{append} is amortized $O(\log \Sigma)$, so
 * building is $O(n \log \Sigma)$. \texttt{order}, \texttt{occ},
 * \texttt{distinct} are $O(n)$, \texttt{kth} is
 * $O(n + |ans| \Sigma)$, \texttt{lcs} is $O(|t| \log \Sigma)$.
 * Usage:
 *  SA sa; int p = 0;
 *  for (char c : s) p = sa.append(p, c);
 *  // for another string, start again from p = 0
 *  vi cnt = sa.occ(s); // cnt[p] = 1 here (p = whole s)
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
	vi order() {
		vi c(N + 1), o(N);
		rep(i,0,N) c[dis[i] + 1]++;
		rep(i,0,N) c[i + 1] += c[i];
		rep(i,0,N) o[c[dis[i]]++] = i;
		return o;
	}
	vi occ(const string& s) {
		vi r(N), o = order();
		int p = 0;
		for (char c : s) r[p = adj[p][c]]++;
		for (int i = N; --i;) r[link[o[i]]] += r[o[i]];
		return r;
	}
	ll distinct() {
		ll r = 0;
		rep(i,1,N) r += dis[i] - dis[link[i]];
		return r;
	}
	string kth(ll k) {
		vi o = order();
		vector<ll> d(N, 1);
		for (int i = N; i--;)
			for (auto [c, v] : adj[o[i]]) d[o[i]] += d[v];
		string r;
		int p = 0;
		for (k++; k--;) for (auto [c, v] : adj[p]) {
			if (k < d[v]) { r += c, p = v; break; }
			k -= d[v];
		}
		return r;
	}
	pii lcs(const string& t) {
		pii r;
		int p = 0, l = 0;
		rep(i,0,sz(t)) {
			while (p && !adj[p].count(t[i])) l = dis[p = link[p]];
			if (adj[p].count(t[i])) p = adj[p][t[i]], l++;
			r = max(r, pii(l, i + 1 - l));
		}
		return r;
	}
};
