#include "../utilities/template.h"
#define pb push_back

#include "../../content/number-theory/ModInt.h"
#include "../../content/numerical/PolyOps.h"
#include "../../content/numerical/PolyMultipoint.h"

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
poly naiveMul(const poly& a, const poly& b, int n) {
	poly c(n);
	rep(i,0,min(n,sz(a))) rep(j,0,min(sz(b),n-i)) c[i+j] += a[i]*b[j];
	return c;
}

int main() {
	rep(it,0,1000) {
		poly a = rp(rnd(0, 30)), b = rp(rnd(0, 30)), c;
		if (sz(a) && sz(b)) c.resize(sz(a)+sz(b)-1);
		rep(i,0,sz(a)) rep(j,0,sz(b)) c[i+j] += a[i]*b[j];
		assert(eq(a*b, c));
	}
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
	rep(it,0,3000) {
		int n = rnd(0, 40); poly a = rp(rnd(1, 45)), one(n);
		if (n) one[0] = 1;
		if (!a[0].v) a[0] = 1;
		assert(eq(naiveMul(a, inv(a, n), n), one));
		a[0] = 1;
		assert(eq(exp(log(a, n), n), RSZ(a, n)));
		a[0] = 0;
		poly e = exp(a, n);
		assert(eq(RSZ(dif(e), max(n-1, 0)), naiveMul(dif(RSZ(a, n)), e, max(n-1, 0))));
		assert(eq(log(e, n), RSZ(a, n)));
	}
	rep(it,0,5000) { // sqrt, with leading zeros and non-residues
		int n = rnd(0, 30); poly a = rp(rnd(0, 20));
		rep(i,0,min(rnd(0, 4), sz(a))) a[i] = 0;
		poly sq = naiveMul(a, a, 2*sz(a)+1), r = sqrt(sq, n);
		assert(sz(r) == n && eq(naiveMul(r, r, n), RSZ(sq, n)));
		// a random series has a root iff its first nonzero term is an even power with a residue coefficient
		a.resize(n); r = sqrt(a, n);
		int t = 0; while (t < n && !a[t].v) t++;
		bool has = t == n || (t % 2 == 0 && pow(a[t], (MOD-1)/2).v == 1);
		assert(has ? sz(r) == n && eq(naiveMul(r, r, n), a) : r.empty());
	}
	rep(it,0,3000) {
		poly f = rp(rnd(0, 40)), g = rp(rnd(1, 40));
		if (!g.back().v) g.back() = 1;
		auto [q, r] = quoRem(f, g);
		assert(sz(r) < sz(g) && sz(q) == max(sz(f)-sz(g)+1, 0));
		assert(eq(REMZ(naiveMul(q, g, sz(f)+1)+r), REMZ(f)));
	}
	rep(it,0,2000) {
		poly f = rp(rnd(0, 40)), p = rp(rnd(0, 40)), r = multiEval(f, p);
		assert(sz(r) == sz(p));
		rep(i,0,sz(p)) assert(r[i].v == eval(f, p[i]).v);
		set<int> xs; while (sz(xs) < sz(p)) xs.insert(rnd(0, MOD-1));
		vector<pair<T,T>> v; for (int x : xs) v.pb({x, rnd(0, MOD-1)});
		shuffle(all(v), rng);
		poly g = interpolate(v);
		assert(sz(g) <= sz(v));
		for (auto [x, y] : v) assert(eval(g, x).v == y.v);
	}
	rep(it,0,3000) { // berlekampMassey recovers a random recurrence
		int d = rnd(0, 8), n = 2*d + rnd(0, 5);
		poly c = rp(d+1), s = rp(d); s.resize(n+30);
		rep(i,d,n+30) rep(j,1,d+1) s[i] += c[j]*s[i-j];
		poly got = berlekampMassey(RSZ(s, n));
		int L = sz(got)-1;
		assert(L <= d && got[0].v == 0);
		rep(i,L,n+30) {
			T x = 0; rep(j,1,L+1) x += got[j]*s[i-j];
			assert(x.v == s[i].v);
		}
		int k = rnd(0, n+29);
		assert(solve_linrec(s, got, L, k).v == s[k].v);
	}
	cout << "Tests passed!" << endl;
}
