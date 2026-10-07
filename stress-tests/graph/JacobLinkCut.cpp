#include "../utilities/template.h"

#include "../../content/graph/JacobLinkCut.h"

mt19937 rng(12345);
int ri(int a, int b) { return (int)(rng() % (unsigned)(b - a + 1)) + a; }

// Naive forest: adjacency sets, path found by DFS.
struct Naive {
	int n; vector<set<int>> adj; vector<ll> val;
	Naive(int n) : n(n), adj(n), val(n) {}
	bool path(int u, int t, int p, vi& out) {
		out.push_back(u);
		if (u == t) return 1;
		for (int v : adj[u]) if (v != p && path(v, t, u, out)) return 1;
		out.pop_back();
		return 0;
	}
	bool conn(int u, int v) { vi p; return path(u, v, -1, p); }
	ll sum(int u, int v) {
		vi p; assert(path(u, v, -1, p));
		ll s = 0; for (int x : p) s += val[x];
		return s;
	}
};

// shape: 0 random, 1 path-ish (link i,i+1), 2 star-ish
void test(int maxn, int iters, int ops, ll maxv) {
	rep(it,0,iters) {
		int n = ri(1, maxn), shape = ri(0, 3);
		vector<node> lct(n);
		Naive nv(n);
		int root = -1; // root of the tree containing `rootOf`, if known
		rep(op,0,ops) {
			int t = ri(0, 9), a = ri(0, n-1), b = ri(0, n-1);
			if (shape == 1 && a + 1 < n && ri(0, 2)) b = a + 1;
			if (shape == 2 && ri(0, 2)) b = 0;
			if (t <= 2) { // link
				if (nv.conn(a, b)) continue;
				link(&lct[a], &lct[b]);
				nv.adj[a].insert(b), nv.adj[b].insert(a);
			} else if (t == 3) { // cut a random existing edge
				if (nv.adj[a].empty()) continue;
				auto e = nv.adj[a].begin();
				advance(e, ri(0, sz(nv.adj[a]) - 1));
				b = *e;
				if (ri(0, 1)) swap(a, b);
				cut(&lct[a], &lct[b]);
				nv.adj[a].erase(b), nv.adj[b].erase(a);
			} else if (t <= 5) { // point add (possibly negative)
				ll v = (ll)(rng() % (unsigned long long)(2 * maxv + 1)) - maxv;
				update(&lct[a], v);
				nv.val[a] += v;
			} else if (t <= 7) { // path sum
				if (!nv.conn(a, b)) continue;
				assert(query(&lct[a], &lct[b]) == nv.sum(a, b));
			} else if (t == 8) { // connectivity via find_root
				bool same = find_root(&lct[a]) == find_root(&lct[b]);
				assert(same == nv.conn(a, b));
			} else { // evert makes a the root of its tree
				evert(&lct[a]); root = a;
				rep(i,0,n) if (nv.conn(a, i))
					assert(find_root(&lct[i]) == &lct[root]);
			}
		}
		// final full check of all pairs
		rep(a,0,n) rep(b,0,n) {
			bool c = nv.conn(a, b);
			assert((find_root(&lct[a]) == find_root(&lct[b])) == c);
			if (c) assert(query(&lct[a], &lct[b]) == nv.sum(a, b));
		}
	}
}

int main() {
	test(1, 200, 20, 10);
	test(2, 2000, 30, 10);
	test(5, 20000, 60, 10);
	test(12, 5000, 200, 1000000000);
	test(40, 300, 1000, (ll)1e13);
	cout<<"Tests passed!"<<endl;
}
