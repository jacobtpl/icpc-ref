#include "../utilities/template.h"

#include "../../content/geometry/Point.h"
#include "../../content/geometry/ManhattanMST.h"
#include "../../content/data-structures/UnionFind.h"


typedef Point<int> P;
typedef int T;
T rectilinear_mst_n(vector<P> ps) {
  struct edge { int src, dst; T weight; };
  vector<edge> edges;

  auto dist = [&](int i, int j) {
    return abs((ps[i]-ps[j]).x) + abs((ps[i]-ps[j]).y);
  };
  for (int i = 0; i < sz(ps); ++i)
    for (int j = i+1; j < sz(ps); ++j)
      edges.push_back({i, j, dist(i,j)});
  T cost = 0;
  sort(all(edges), [](edge a, edge b) { return a.weight < b.weight; });
  UF uf(sz(ps));
  for (auto e: edges)
    if (uf.join(e.src, e.dst))
      cost += e.weight;
  return cost;
}

// O(N^2) Prim with 64-bit weights
ll primCost(const vector<P>& ps) {
	int n = sz(ps);
	if (!n) return 0;
	vector<ll> d(n, LLONG_MAX); vi done(n);
	ll cost = 0; d[0] = 0;
	rep(it,0,n) {
		int b = -1;
		rep(i,0,n) if (!done[i] && (b == -1 || d[i] < d[b])) b = i;
		done[b] = 1; cost += d[b];
		rep(i,0,n) d[i] = min(d[i],
			abs((ll)ps[i].x - ps[b].x) + abs((ll)ps[i].y - ps[b].y));
	}
	return cost;
}

void check(const vector<P>& pts) {
	auto edges = manhattanMST(pts);
	assert(sz(edges) <= 4*sz(pts));
	for (auto e : edges) {
		assert(0 <= e[1] && e[1] < sz(pts) && 0 <= e[2] && e[2] < sz(pts));
		P d = pts[e[1]] - pts[e[2]];
		assert(e[0] == abs(d.x) + abs(d.y));
	}
	sort(all(edges));
	UF uf(sz(pts));
	ll cost = 0; int joined = 0;
	for (auto e: edges) if (uf.join(e[1], e[2])) cost += e[0], joined++;
	assert(joined == max(sz(pts) - 1, 0));
	assert(cost == primCost(pts));
}

signed main() {
	srand(5);
	check({});
	check({P(3,-4)});
	check({P(3,-4), P(3,-4)});
	check({P(0,0), P(1000000000,1000000000)});
	check({P(-500000000,-500000000), P(500000000,500000000),
			P(-500000000,500000000), P(500000000,-500000000), P(0,0)});
#ifdef TEST_INT_OVERFLOW
	// Distances that do not fit in an int: the returned weight wraps around
	// (4e9 becomes -294967296). Opt-in, since this needs P/edges to be ll.
	check({P(-1000000000,-1000000000), P(1000000000,1000000000)});
#endif
	// tiny grids: duplicates, ties, collinear points
	rep(t,0,300000) {
		int n = rand() % 9, C = rand() % 4 + 1;
		vector<P> pts;
		rep(i,0,n) pts.push_back(P(rand() % C - C / 2, rand() % C - C / 2));
		check(pts);
	}
	// structured and large-coordinate inputs. All pairwise distances must fit
	// in an int, so coordinates stay within [-5e8, 5e8].
	rep(t,0,3000) {
		int n = rand() % 200, type = t % 7, C = t % 2 ? 500000000 : 50;
		vector<P> pts;
		rep(i,0,n) {
			int x = rand() % (2*C+1) - C, y = rand() % (2*C+1) - C;
			if (type == 1) y = x;  // main diagonal
			if (type == 2) y = -x; // anti-diagonal
			if (type == 3) y = 7;  // horizontal line
			if (type == 4) x = -3; // vertical line
			if (type == 5) x = x / (C/5) * (C/5), y = y / (C/5) * (C/5); // grid
			if (type == 6) // clusters in the four corners
				x = (x > 0 ? C : -C) - x % 3, y = (y > 0 ? C : -C) - y % 3;
			pts.push_back(P(x, y));
		}
		check(pts);
	}
    for (int t=0; t<10000; t++) {
        const int max_coord = rand() % 300 + 1;
        const int num_pts = rand() % 100;
        vector<P> pts;
        for (int i = 0; i < num_pts; ++i) {
            int x = rand() % max_coord - max_coord / 2;
            int y = rand() % max_coord - max_coord / 2;
            pts.push_back(P(x,y));
        }
        auto edges = manhattanMST(pts);
        assert(edges.size() <= 4*pts.size());
        sort(all(edges));
        UF uf(sz(pts));
        int cost = 0, joined = 0;
        for (auto e: edges) if (uf.join(e[1], e[2])) cost += e[0], joined++;
        if (num_pts > 0) assert(joined == num_pts - 1);
        assert(cost == rectilinear_mst_n(pts));
    }
    cout<<"Tests passed!"<<endl;
}
