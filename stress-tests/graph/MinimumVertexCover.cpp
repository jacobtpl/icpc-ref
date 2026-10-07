#include "../utilities/template.h"

#include "../../content/graph/MinimumVertexCover.h"
#include "../../content/graph/hopcroftKarp.h"

vi coverHK(vector<vi>& g, int n, int m) {
	vi match(m, -1);
	int res = hopcroftKarp(g, match);
	vector<bool> lfound(n, true), seen(m);
	for(auto &it: match) if (it != -1) lfound[it] = false;
	vi q, cover;
	rep(i,0,n) if (lfound[i]) q.push_back(i);
	while (!q.empty()) {
		int i = q.back(); q.pop_back();
		lfound[i] = 1;
		for(auto &e: g[i]) if (!seen[e] && match[e] != -1) {
			seen[e] = true;
			q.push_back(match[e]);
		}
	}
	rep(i,0,n) if (!lfound[i]) cover.push_back(i);
	rep(i,0,m) if (seen[i]) cover.push_back(n+i);
	assert(sz(cover) == res);
	return cover;
}

int main() {
	rep(it,0,300000) {
		int N = rand() % 20, M = rand() % 20;
		int prop = rand();
		vector<vi> gr(N);
		vi left(N), right(M);
		rep(i,0,N) rep(j,0,M) if (rand() < prop) {
			gr[i].push_back(j);
		}
		auto verify = [&](vi& cover) {
			for(auto &x: cover) {
				if (x < N) left[x] = 1;
				else right[x - N] = 1;
			}
			rep(i,0,N) if (!left[i]) for(auto &j:gr[i]) {
				assert(right[j]);
				/* if (!right[j]) {
					cout << N << ' ' << M << endl;
					rep(i,0,N) for(auto &j: gr[i]) cout << i << " - " << j << endl;
					cout << "yields " << sz(cover) << endl;
					for(auto &x: cover) cout << x << endl;
					abort();
				} */
			}
		};
		vi cover1 = cover(gr, N, M);
		vi cover2 = coverHK(gr, N, M);
		assert(sz(cover1) == sz(cover2));
		verify(cover1);
		verify(cover2);
		// cout << '.' << endl;
	}
	// exact minimum by brute force over all vertex subsets (multi-edges, empty sides)
	{ vector<vi> g; assert(cover(g, 0, 0).empty()); assert(cover(g, 0, 3).empty()); }
	{ vector<vi> g(3); assert(cover(g, 3, 0).empty()); }
	mt19937 rng(99);
	rep(it,0,30000) {
		int N = (int)(rng() % 7), M = (int)(rng() % 7), p = (int)(rng() % 101);
		vector<vi> gr(N);
		rep(i,0,N) rep(j,0,M) rep(k,0,2) if ((int)(rng() % 100) < (k ? p / 3 : p)) gr[i].push_back(j);
		rep(i,0,N) shuffle(all(gr[i]), rng);
		int best = N + M;
		rep(mask,0,1 << (N + M)) {
			bool ok = 1;
			rep(i,0,N) if (!(mask >> i & 1)) for (int j : gr[i]) if (!(mask >> (N + j) & 1)) ok = 0;
			if (ok) best = min(best, __builtin_popcount(mask));
		}
		for (vi c : {cover(gr, N, M), coverHK(gr, N, M)}) {
			assert(sz(c) == best);
			int mask = 0;
			for (int x : c) { assert(0 <= x && x < N + M && !(mask >> x & 1)); mask |= 1 << x; }
			rep(i,0,N) if (!(mask >> i & 1)) for (int j : gr[i]) assert(mask >> (N + j) & 1);
		}
	}
	cout<<"Tests passed!"<<endl;
}
