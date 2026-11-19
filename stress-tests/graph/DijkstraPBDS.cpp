#include "../utilities/template.h"

#include "../../content/graph/DijkstraPBDS.h"
#include "../../content/graph/BellmanFord.h"

struct HopsDist {
	typedef pair<ll, int> D;
	D zero() { return {0, 0}; }
	D add(D d, ll w) { return {d.first + w, d.second + 1}; }
	bool less(D a, D b) { return a < b; }
};

template<class D> void checkPar(const DijkstraRes<D>& r, const vector<vector<pair<int, ll>>>& g,
		vi src, const vector<ll>& want) {
	rep(i,0,sz(g)) if (r.seen[i]) {
		int v = i; ll len = 0; int steps = 0;
		while (r.par[v] != -1) {
			int p = r.par[v]; ll best = LLONG_MAX;
			for (auto [x, w] : g[p]) if (x == v) best = min(best, w);
			assert(best != LLONG_MAX);
			len += best; v = p; assert(++steps <= sz(g));
		}
		assert(count(all(src), v) && len == want[i]);
	}
}

int main() {
	rep(it,0,20000) {
		int n = rand() % (it < 10000 ? 8 : 60) + 1, m = rand() % (4 * n);
		int W = rand() % 3 ? 10 : 1000000;
		vector<vector<pair<int, ll>>> g(n);
		vector<vector<pii>> gi(n);
		vector<Ed> eds;
		rep(i,0,m) {
			int a = rand() % n, b = rand() % n, w = rand() % W;
			g[a].emplace_back(b, w); gi[a].emplace_back(b, w);
			eds.push_back({a, b, w});
		}
		int s = rand() % n, s2 = rand() % n;
		vector<Node> n1(n), n2(n);
		bellmanFord(n1, eds, s); bellmanFord(n2, eds, s2);
		vector<ll> want(n), want2(n);
		rep(i,0,n) want[i] = n1[i].dist, want2[i] = min(n1[i].dist, n2[i].dist);

		auto r0 = dijkstraPBDS(g, {s});
		auto r1 = dijkstraPBDS(gi, {s});
		auto rh = dijkstraPBDS(g, {s}, HopsDist{});
		auto rd = dijkstra(g, {s}, HopsDist{});
		auto r2 = dijkstraPBDS(g, {s, s2});
		auto r3 = dijkstraPBDSGen(n, {s}, LLDist{}, [&](int u, auto f) {
			for (auto [v, w] : g[u]) f(v, w);
		});
		rep(i,0,n) {
			bool reach = want[i] != inf;
			assert(r0.seen[i] == reach && r1.seen[i] == reach && rh.seen[i] == reach && r3.seen[i] == reach);
			assert((bool)r2.seen[i] == (want2[i] != inf));
			if (reach) assert(r0.dist[i] == want[i] && r1.dist[i] == want[i] && r3.dist[i] == want[i]);
			if (reach) assert(rh.dist[i] == rd.dist[i]); // same (length, hops) as Dijkstra.h
			if (r2.seen[i]) assert(r2.dist[i] == want2[i]);
		}
		checkPar(r0, g, {s}, want);
		checkPar(r2, g, {s, s2}, want2);
		rep(i,0,n) if (rh.seen[i]) {
			int v = i, hops = 0;
			while (rh.par[v] != -1) v = rh.par[v], hops++;
			assert(hops == rh.dist[i].second);
		}
	}
	// large: every vertex is decreased many times (heap modify path)
	{
		int n = 2000;
		vector<vector<pair<int, ll>>> g(n);
		ll C = (ll)n * n + 2 * n;
		rep(i,0,n-1) {
			g[i].emplace_back(i + 1, 1);
			rep(j,i+2,n) g[i].emplace_back(j, C - (ll)i * n - i + rand() % n);
		}
		auto a = dijkstraPBDS(g, {0});
		auto b = dijkstra(g, {0});
		assert(a.dist == b.dist);
		rep(i,0,n) assert(a.dist[i] == i);
	}
	{
		int n = 200000;
		vector<vector<pair<int, ll>>> g(n);
		rep(i,0,n-1) g[i].emplace_back(i + 1, 1000000000);
		auto r = dijkstraPBDS(g, {0});
		assert(r.dist[n-1] == (ll)(n - 1) * 1000000000);
	}
	cout << "Tests passed!" << endl;
}
