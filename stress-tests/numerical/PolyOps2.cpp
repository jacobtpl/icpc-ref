#include "../utilities/template.h"
#define pb push_back

#include "../../content/number-theory/SiyongModular.h"
#include "../../content/numerical/PolyOps.h"
#include "../../content/numerical/PolyOps2.h"

mt19937 rng(2);
int rnd(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }
poly rp(int n) { poly p(n); for (auto& x : p) x = rnd(0, MOD-1); return p; }
bool eq(const poly& a, const poly& b) {
	if (sz(a) != sz(b)) return 0;
	rep(i,0,sz(a)) if (a[i].v != b[i].v) return 0;
	return 1;
}
poly naivePow(poly a, ll b, int n) {
	poly r(n); r[0] = 1; a.resize(n);
	for (; b; b /= 2, a = RSZ(a*a, n)) if (b & 1) r = RSZ(r*a, n);
	return r;
}

int main() {
	rep(it,0,5000) {
		int n = rnd(1, 12);
		poly a = rp(rnd(0, n+2));
		rep(i,0,min(rnd(0, 3), sz(a))) a[i] = 0;
		ll b = rnd(0, 2) ? rnd(0, 8) : (ll)rnd(0, 1e9) * rnd(0, 1e9);
		assert(eq(pow(a, b, n), naivePow(a, b, n)));
	}
	rep(it,0,1000) {
		int d = rnd(1, 8); ll k = rnd(0, 2) ? rnd(0, 40) : (ll)rnd(0, 1e9) * rnd(0, 1e9);
		poly s = rp(d), c = rp(d+1);
		T res = solve_linrec(s, c, d, k);
		if (k <= 40) {
			poly seq = s; seq.resize(k+d+1);
			rep(i,d,k+1) rep(j,1,d+1) seq[i] += c[j]*seq[i-j];
			assert(res.v == seq[k].v);
		}
		poly f(d+1); f[d] = 1; rep(i,0,d) f[i] = -c[d-i];
		poly r = xkmodf(k, f); r.resize(d);
		T ans = 0; rep(i,0,d) ans += r[i]*s[i];
		assert(res.v == ans.v);

		poly Q = rp(d+1); Q[0] = rnd(1, MOD-1); poly P = rp(rnd(0, d));
		k %= 50;
		assert(kthCoef(P, Q, k).v == RSZ(conv(P, inv(Q, int(k)+1)), int(k)+1)[k].v);
	}
	rep(it,0,1000) {
		poly f = rp(rnd(0, 20)); T c = rnd(0, MOD-1), x = rnd(0, MOD-1);
		poly g = taylorShift(f, c);
		assert(sz(g) == sz(f) && eval(g, x).v == eval(f, x+c).v);
	}
	cout << "Tests passed!" << endl;
}
