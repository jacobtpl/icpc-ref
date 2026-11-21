#include "../utilities/template.h"

#include "../../content/data-structures/OfflineDynamicConnectivity.h"

// Brute force: multiset of edges, DFS at every query.
struct Brute {
	int n;
	multiset<pii> es;
	Brute(int n) : n(n) {}
	void add(int a, int b) { es.insert(minmax(a, b)); }
	void rem(int a, int b) {
		auto it = es.find(minmax(a, b));
		assert(it != es.end());
		es.erase(it);
	}
	vi comps() {
		vector<vi> adj(n);
		for (auto& e : es) {
			adj[e.first].push_back(e.second);
			adj[e.second].push_back(e.first);
		}
		vi c(n, -1);
		rep(i,0,n) if (c[i] < 0) {
			vi st = {i}; c[i] = i;
			while (!st.empty()) {
				int x = st.back(); st.pop_back();
				for (int y : adj[x]) if (c[y] < 0) c[y] = i, st.push_back(y);
			}
		}
		return c;
	}
};

void check(DynCon& d, vector<vi>& expected, int n) {
	int calls = 0;
	d.solve([&](int i) {
		assert(i == calls++);
		vi& c = expected[i];
		int ncomp = 0;
		rep(x,0,n) {
			if (c[x] == x) ncomp++;
			int s = 0;
			rep(y,0,n) {
				assert((c[x] == c[y]) == (d.uf.find(x) == d.uf.find(y)));
				if (c[x] == c[y]) s++;
			}
			assert(d.uf.size(x) == s);
		}
		assert(ncomp == n - sz(d.uf.st) / 2);
	});
	assert(calls == sz(expected));
	assert(d.uf.st.empty());
	rep(x,0,n) assert(d.uf.e[x] == -1);
}

void randomTest(int maxN, int maxOps, int iters) {
	rep(it,0,iters) {
		int n = rand() % maxN + 1;
		int ops = rand() % (maxOps + 1);
		int mode = rand() % 4; // 0: any pair, 1: few pairs (multi-edges), 2: never remove, 3: remove-heavy
		int pairs = mode == 1 ? rand() % 3 + 1 : 0;
		DynCon d(n);
		Brute b(n);
		vector<vi> expected;
		vector<pii> alive;
		rep(o,0,ops) {
			int t = rand() % 10;
			if (t < 3 || (mode == 2 && t < 7)) {
				int a, c;
				if (pairs) { int p = rand() % pairs; a = p % n; c = (p * 7 + 3) % n; }
				else { a = rand() % n; c = rand() % n; }
				d.add(a, c); b.add(a, c); alive.push_back({a, c});
			} else if (t < (mode == 3 ? 8 : 6) && !alive.empty() && mode != 2) {
				int k = rand() % sz(alive);
				pii e = alive[k];
				alive.erase(alive.begin() + k);
				if (rand() % 2) swap(e.first, e.second);
				d.rem(e.first, e.second); b.rem(e.first, e.second);
			} else {
				int i = d.query();
				assert(i == sz(expected));
				expected.push_back(b.comps());
			}
		}
		check(d, expected, n);
	}
}

void edgeCases() {
	{ // n = 1, no queries
		DynCon d(1);
		vector<vi> e;
		check(d, e, 1);
	}
	{ // n = 1, self-loops, queries
		DynCon d(1);
		Brute b(1);
		vector<vi> e;
		d.add(0, 0); b.add(0, 0); d.query(); e.push_back(b.comps());
		d.rem(0, 0); b.rem(0, 0); d.query(); e.push_back(b.comps());
		check(d, e, 1);
	}
	{ // no queries, only edges
		DynCon d(5);
		vector<vi> e;
		d.add(0, 1); d.add(1, 2); d.rem(2, 1);
		check(d, e, 5);
	}
	{ // edge added and removed with no query in between; edge added after last query
		DynCon d(3);
		Brute b(3);
		vector<vi> e;
		d.query(); e.push_back(b.comps());
		d.add(0, 1); d.rem(1, 0);
		d.query(); e.push_back(b.comps());
		d.add(1, 2); b.add(1, 2);
		d.query(); e.push_back(b.comps());
		d.add(0, 2);
		check(d, e, 3);
	}
	{ // duplicate edge removed one copy at a time
		DynCon d(2);
		Brute b(2);
		vector<vi> e;
		d.add(0, 1); b.add(0, 1); d.add(1, 0); b.add(1, 0);
		d.query(); e.push_back(b.comps());
		d.rem(0, 1); b.rem(0, 1);
		d.query(); e.push_back(b.comps());
		d.rem(1, 0); b.rem(1, 0);
		d.query(); e.push_back(b.comps());
		check(d, e, 2);
	}
	{ // large-ish: path built then torn down, many queries
		int n = 300;
		DynCon d(n);
		Brute b(n);
		vector<vi> e;
		rep(i,0,n-1) { d.add(i, i + 1); b.add(i, i + 1); if (i % 10 == 0) { d.query(); e.push_back(b.comps()); } }
		rep(i,0,n-1) { d.rem(i + 1, i); b.rem(i + 1, i); if (i % 10 == 0) { d.query(); e.push_back(b.comps()); } }
		check(d, e, n);
	}
}

int main() {
	edgeCases();
	randomTest(7, 40, 30000);
	randomTest(3, 60, 10000);
	randomTest(40, 400, 300);
	cout<<"Tests passed!"<<endl;
}
