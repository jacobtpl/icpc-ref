#include "../utilities/template.h"

#include "../../content/data-structures/UnionFind.h"

mt19937 rng(4242);
int rnd(int a, int b) { return a + (int)(rng() % (unsigned)(b - a + 1)); }

void run(int n, int ops) {
	UF uf(n);
	vi comp(n);
	iota(all(comp), 0);
	rep(it,0,ops) {
		int a = rnd(0, n - 1), b = rnd(0, 3) ? rnd(0, n - 1) : a; // self-joins too
		int t = rnd(0, 3);
		if (t == 0) {
			bool exp = comp[a] != comp[b];
			int ca = comp[a], cb = comp[b];
			rep(i,0,n) if (comp[i] == cb) comp[i] = ca;
			assert(uf.join(a, b) == exp);
		} else if (t == 1) {
			assert(uf.sameSet(a, b) == (comp[a] == comp[b]));
		} else if (t == 2) {
			assert(uf.size(a) == (int)count(all(comp), comp[a]));
		} else {
			int f = uf.find(a);
			assert(0 <= f && f < n && comp[f] == comp[a]);
			assert(uf.find(f) == f && (uf.find(b) == f) == (comp[a] == comp[b]));
		}
	}
	rep(i,0,n) rep(j,0,n) assert(uf.sameSet(i, j) == (comp[i] == comp[j]));
}

int main() {
	{ UF uf(0); UF u1(1); assert(u1.find(0) == 0 && u1.size(0) == 1 && !u1.join(0, 0)); }
	rep(it,0,30000) run(rnd(1, 8), rnd(0, 30));
	rep(it,0,300) run(rnd(1, 200), rnd(0, 600));
	{ // chains built in both directions stay shallow and sizes add up
		int n = 1 << 20;
		UF uf(n);
		rep(i,1,n) assert(uf.join(i - 1, i));
		assert(uf.size(0) == n && uf.sameSet(0, n - 1));
		UF u2(n);
		for (int i = n - 1; i > 0; i--) assert(u2.join(i, i - 1));
		assert(u2.size(n - 1) == n && !u2.join(0, n - 1));
	}
	cout<<"Tests passed!"<<endl;
}
