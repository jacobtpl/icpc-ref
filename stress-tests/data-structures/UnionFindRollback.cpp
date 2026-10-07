#include "../utilities/template.h"

#include "../../content/data-structures/UnionFindRollback.h"

mt19937 rng(9001);
int rnd(int a, int b) { return a + (int)(rng() % (unsigned)(b - a + 1)); }

vi components(int n, const vector<pii>& ed) {
	vi comp(n);
	iota(all(comp), 0);
	for (auto [a, b] : ed) {
		int ca = comp[a], cb = comp[b];
		rep(i,0,n) if (comp[i] == cb) comp[i] = ca;
	}
	return comp;
}

void run(int n, int ops) {
	RollbackUF uf(n);
	vector<pii> ed; // successful joins, in order
	vi comp = components(n, ed);
	vi times = {0}; // times[i] = uf.time() after i successful joins
	rep(it,0,ops) {
		int a = rnd(0, n - 1), b = rnd(0, 3) ? rnd(0, n - 1) : a;
		int t = rnd(0, 5);
		if (t <= 1) {
			bool exp = comp[a] != comp[b];
			assert(uf.join(a, b) == exp);
			if (exp) {
				ed.push_back({a, b}), times.push_back(uf.time());
				assert(times[sz(times) - 2] < times.back());
			}
			assert(uf.time() == times.back());
			comp = components(n, ed);
		} else if (t == 2) {
			int k = rnd(0, sz(ed)); // also a no-op rollback to the current time
			uf.rollback(times[k]);
			ed.resize(k), times.resize(k + 1);
			assert(uf.time() == times.back());
			comp = components(n, ed);
		} else if (t == 3) {
			assert(uf.size(a) == (int)count(all(comp), comp[a]));
		} else {
			int f = uf.find(a);
			assert(0 <= f && f < n && comp[f] == comp[a]);
			assert((uf.find(b) == f) == (comp[a] == comp[b]));
		}
	}
	rep(i,0,n) rep(j,0,n) assert((uf.find(i) == uf.find(j)) == (comp[i] == comp[j]));
	uf.rollback(0);
	rep(i,0,n) assert(uf.find(i) == i && uf.size(i) == 1);
}

int main() {
	{ RollbackUF uf(0); uf.rollback(0); RollbackUF u1(1); assert(!u1.join(0, 0) && u1.time() == 0); }
	rep(it,0,30000) run(rnd(1, 8), rnd(0, 40));
	rep(it,0,300) run(rnd(1, 100), rnd(0, 500));
	{ // find depth stays logarithmic on chains, in both directions
		int n = 1 << 20;
		RollbackUF uf(n);
		rep(i,1,n) assert(uf.join(i - 1, i));
		assert(uf.size(0) == n);
		uf.rollback(0);
		for (int i = n - 1; i > 0; i--) assert(uf.join(i, i - 1));
		assert(uf.size(n - 1) == n && !uf.join(0, n - 1));
	}
	cout<<"Tests passed!"<<endl;
}
