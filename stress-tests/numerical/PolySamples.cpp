#include "../utilities/template.h"

#include "../../content/number-theory/ModInt.h"
#include "../../content/numerical/PolyOps.h"
#include "../../content/numerical/PolySamples.h"

mt19937 rng(3);
int rnd(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }

int main() {
	rep(it,0,20000) {
		int n = rnd(0, 12), m = rnd(0, 15);
		poly f(n), y(n);
		for (auto& x : f) x = rnd(0, MOD-1);
		rep(i,0,n) y[i] = eval(f, i);
		// c near 0 and near MOD so that the points overlap 0..n-1
		T c = rnd(0, 2) ? rnd(-20, 20) : rnd(0, MOD-1);
		poly r = shiftSamples(y, c, m);
		assert(sz(r) == m);
		rep(j,0,m) assert(r[j].v == eval(f, c+j).v);
		assert(lagrangeAt(y, c).v == eval(f, c).v);
	}
	cout << "Tests passed!" << endl;
}
