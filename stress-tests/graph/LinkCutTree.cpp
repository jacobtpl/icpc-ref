#include "../utilities/template.h"

#include "../../content/graph/LinkCutTree.h"
#include "../../content/data-structures/UnionFind.h"

int main() {
	srand(2);
	LinkCut lczero(0);
	rep(it,0,10000) {
		int N = rand() % 20 + 1;
		LinkCut lc(N);
		UF uf(N);
		vector<pii> edges;
		rep(it2,0,1000) {
			int v = (rand() >> 4) & 3;
			if (v == 0 && !edges.empty()) { // remove
				int r = (rand() >> 4) % sz(edges);
				pii ed = edges[r];
				swap(edges[r], edges.back());
				edges.pop_back();
				if (rand() & 16)
					lc.cut(ed.first, ed.second);
				else
					lc.cut(ed.second, ed.first);
			} else {
				int a = (rand() >> 4) % N;
				int b = (rand() >> 4) % N;
				uf.e.assign(N, -1);
				for(auto &ed: edges) uf.join(ed.first, ed.second);
				bool c = uf.sameSet(a, b);
				if (!c && v != 1) {
					lc.link(a, b);
					edges.emplace_back(a, b);
				} else {
					assert(lc.connected(a, b) == c);
				}
			}
		}
	}
	// Larger forests, including long paths and stars; oracle is a BFS.
	rep(it,0,60) {
		int N = rand() % 300 + 2, shape = it % 3;
		LinkCut lc(N);
		vector<set<int>> adj(N);
		vector<pii> edges;
		auto conn = [&](int a, int b) {
			vi seen(N), q = {a}; seen[a] = 1;
			rep(i,0,sz(q)) for (int y : adj[q[i]])
				if (!seen[y]) seen[y] = 1, q.push_back(y);
			return (bool)seen[b];
		};
		auto add = [&](int a, int b) {
			lc.link(a, b); adj[a].insert(b); adj[b].insert(a);
			edges.emplace_back(a, b);
		};
		if (shape == 1) rep(i,1,N) add(i - 1, i);
		if (shape == 2) rep(i,1,N) rand() & 16 ? add(0, i) : add(i, 0);
		rep(it2,0,3000) {
			int a = (rand() >> 4) % N, b = (rand() >> 4) % N;
			bool c = conn(a, b);
			assert(lc.connected(a, b) == c);
			assert(lc.connected(a, a));
			if (!c) add(a, b);
			else if (!edges.empty() && rand() % 3) {
				int r = (rand() >> 4) % sz(edges);
				tie(a, b) = edges[r];
				edges[r] = edges.back(); edges.pop_back();
				adj[a].erase(b); adj[b].erase(a);
				if (rand() & 16) swap(a, b);
				lc.cut(a, b);
				assert(!lc.connected(a, b));
			}
		}
	}
	cout<<"Tests passed!"<<endl;
}
