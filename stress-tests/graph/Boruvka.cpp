#include "../utilities/template.h"

#include "../../content/graph/Boruvka.h"

int main() {
	rep(it,0,50000) {
		int n = rand() % 9 + 1, m = rand() % 20, W = rand() % 2 ? 3 : 100;
		vector<tuple<ll,int,int>> es;
		rep(i,0,m) {
			int u = rand() % n, v = rand() % n;
			if (u != v) es.emplace_back(rand() % W - W / 2, u, v);
		}
		int calls = 0;
		auto res = boruvka<ll>(n, [&](vi& comp, vector<tuple<ll,int,int>>& best) {
			calls++;
			for (auto [w, u, v] : es) rep(t,0,2) {
				if (comp[u] != comp[v]) {
					auto &b = best[comp[u]];
					if (get<2>(b) < 0 || w < get<0>(b)) b = {w, u, v};
				}
				swap(u, v);
			}
		});
		sort(all(es));
		UF uf(n); ll tot = 0; int cnt = 0;
		for (auto [w, u, v] : es) if (uf.join(u, v)) tot += w, cnt++;
		assert(res.first == tot && sz(res.second) == cnt);
		UF uf2(n);
		for (auto [u, v] : res.second) assert(uf2.join(u, v));
		int lg = 0; while ((1 << lg) < n) lg++;
		assert(calls <= lg + 1);
	}
	cout<<"Tests passed!"<<endl;
}
