#include "../utilities/template.h"

#include "../../content/graph/BlockCutTree.h"

mt19937 rng(987654);
int ri(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }

int components(int n, const vector<pii>& es, int del, vi& comp) {
	comp.assign(n, -1);
	int c = 0;
	rep(s,0,n) if (comp[s] == -1 && s != del) {
		vi q = {s};
		comp[s] = c;
		rep(i,0,sz(q)) for (auto e : es) rep(k,0,2) {
			if (e.first == q[i] && e.second != del && comp[e.second] == -1)
				comp[e.second] = c, q.push_back(e.second);
			swap(e.first, e.second);
		}
		c++;
	}
	return c;
}

void check(int n, const vector<pii>& es) {
	int m = sz(es);
	ed.assign(n, {});
	edges = es;
	rep(i,0,m) {
		ed[es[i].first].emplace_back(es[i].second, i);
		ed[es[i].second].emplace_back(es[i].first, i);
	}
	for (auto& v : ed) shuffle(all(v), rng);
	auto [cut, tadj, who, emap, vmap] = BCTree();
	int TN = sz(tadj);
	assert(sz(who) == TN && sz(emap) == m && sz(vmap) == n && 0 <= cut && cut <= TN);

	vi comp, tmp;
	int nc = components(n, es, -1, comp);
	vi deg(n);
	for (auto e : es) deg[e.first]++, deg[e.second]++;
	// true articulation points and bridges by brute force
	vector<bool> ap(n), bridge(m);
	rep(v,0,n) ap[v] = components(n, es, v, tmp) > nc - (deg[v] == 0);
	rep(i,0,m) {
		vector<pii> rest = es;
		rest.erase(rest.begin() + i);
		bridge[i] = components(n, rest, -1, tmp) > nc;
	}
	// edges: bridges have emap -1, others a block id below cut
	rep(i,0,m) {
		if (bridge[i]) assert(emap[i] == -1);
		else assert(0 <= emap[i] && emap[i] < cut);
	}
	// two non-bridge edges share a block iff no vertex removal separates them
	rep(i,0,m) rep(j,0,i) if (!bridge[i] && !bridge[j]) {
		bool same = comp[es[i].first] == comp[es[j].first];
		rep(v,0,n) if (same) {
			components(n, es, v, tmp);
			int a = es[i].first != v ? es[i].first : es[i].second;
			int b = es[j].first != v ? es[j].first : es[j].second;
			if (tmp[a] != tmp[b]) same = 0;
		}
		assert(same == (emap[i] == emap[j]));
	}
	// vertices
	vector<bool> bridgeEnd(n);
	rep(i,0,m) if (bridge[i]) bridgeEnd[es[i].first] = bridgeEnd[es[i].second] = 1;
	rep(v,0,n) {
		assert(0 <= vmap[v] && vmap[v] < TN);
		assert(count(all(who[vmap[v]]), v) == 1);
		if (ap[v]) assert(vmap[v] >= cut);
		if (vmap[v] >= cut) {
			// own node: an articulation point, a bridge endpoint or isolated
			assert(ap[v] || bridgeEnd[v] || deg[v] == 0);
			assert(who[vmap[v]] == vi{v});
		} else for (auto [x, e] : ed[v]) assert(emap[e] == vmap[v]);
	}
	int total = 0;
	rep(t,0,TN) {
		total += sz(who[t]);
		assert(t < cut || sz(who[t]) == 1);
	}
	assert(total == n);
	// tree adjacency: simple, symmetric, a forest with one tree per component
	int tedges = 0;
	rep(t,0,TN) {
		assert(is_sorted(all(tadj[t])));
		rep(i,0,sz(tadj[t])) {
			int u = tadj[t][i];
			assert(0 <= u && u < TN && u != t);
			if (i) assert(tadj[t][i-1] != u);
			assert(binary_search(all(tadj[u]), t));
			assert(t >= cut || u >= cut); // blocks are never adjacent
			tedges++;
		}
	}
	vi tcomp(TN, -1);
	int tc = 0;
	rep(s,0,TN) if (tcomp[s] == -1) {
		vi q = {s};
		tcomp[s] = tc;
		rep(i,0,sz(q)) for (int y : tadj[q[i]]) if (tcomp[y] == -1)
			tcomp[y] = tc, q.push_back(y);
		tc++;
	}
	assert(tedges / 2 == TN - tc);
	assert(tc == nc);
	rep(u,0,n) rep(v,0,u)
		assert((comp[u] == comp[v]) == (tcomp[vmap[u]] == tcomp[vmap[v]]));
	// a block node is adjacent exactly to the cut nodes of its vertices
	rep(i,0,m) if (!bridge[i]) for (int v : {es[i].first, es[i].second})
		if (vmap[v] >= cut) assert(binary_search(all(tadj[emap[i]]), vmap[v]));
	rep(i,0,m) if (bridge[i])
		assert(binary_search(all(tadj[vmap[es[i].first]]), vmap[es[i].second]));
}

int main() {
	rep(it,0,60000) {
		int n = ri(1, 8), m = n == 1 ? 0 : ri(0, it % 4 ? 2*n : 3);
		vector<pii> es; // isolated vertices, multi-edges allowed; no self-loops
		rep(i,0,m) {
			int a = ri(0, n-1), b = ri(0, n-2);
			if (b >= a) b++;
			es.emplace_back(a, b);
		}
		check(n, es);
	}
	check(1, {});
	check(3, {});
	check(3, {{0, 1}});
	check(2, {{0, 1}});
	check(2, {{0, 1}, {0, 1}});
	check(3, {{0, 1}, {1, 2}, {2, 0}});
	cout << "Tests passed!" << endl;
}
