// Tests PolyOps.h, PolyOps2.h and PolyMultipoint.h against O(n^2) oracles on plain ll.
#include "../utilities/template.h"
#define pb push_back

#include "../../content/numerical/PolyMultipoint.h"

typedef vector<ll> vl;
const ll M = 998244353;
mt19937_64 rng(2024);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }
ll mpow(ll b, ll e) { ll r = 1; for (b %= M; e; e /= 2, b = b * b % M) if (e & 1) r = r * b % M; return r; }
ll minv(ll a) { return mpow(a, M - 2); }

poly P(const vl& a) { poly p; for (ll x : a) p.pb(mint(x)); return p; }
vl V(const poly& p) { vl a; for (mint x : p) { assert(0 <= x.v && x.v < MOD); a.pb(x.v); } return a; }
ll coef() { int m = (int)rnd(0, 9); return m == 0 ? 0 : m == 1 ? M - 1 : m == 2 ? 1 : rnd(0, M - 1); }
vl randPoly(int n) { vl a(n); for (ll& x : a) x = coef(); return a; }

vl nMul(const vl& a, const vl& b) {
	if (a.empty() || b.empty()) return {};
	vl c(sz(a) + sz(b) - 1);
	rep(i,0,sz(a)) rep(j,0,sz(b)) c[i + j] = (c[i + j] + a[i] * b[j]) % M;
	return c;
}
vl nMulT(const vl& a, const vl& b, int n) { vl c = nMul(a, b); c.resize(n); return c; }
ll nEval(const vl& a, ll x) { ll r = 0; for (int i = sz(a); i--;) r = (r * x + a[i]) % M; return r; }
vl nInv(vl a, int n) {
	a.resize(n); vl b(n); ll i0 = minv(a[0]);
	rep(i,0,n) {
		ll s = i == 0;
		rep(j,1,i+1) s = (s - a[j] * b[i - j]) % M;
		b[i] = (s % M + M) * i0 % M;
	}
	return b;
}
vl nLog(vl a, int n) { // a[0] = 1
	a.resize(n); vl b(n);
	rep(i,1,n) {
		ll s = i * a[i] % M;
		rep(k,1,i) s = (s - k * b[k] % M * a[i - k]) % M;
		b[i] = (s + M) % M * minv(i) % M;
	}
	return b;
}
vl nExp(vl a, int n) { // a[0] = 0
	a.resize(n); vl e(n); if (n) e[0] = 1;
	rep(i,1,n) {
		ll s = 0;
		rep(k,1,i+1) s = (s + k * a[k] % M * e[i - k]) % M;
		e[i] = s * minv(i) % M;
	}
	return e;
}
vl nPow(vl a, ll b, int n) {
	a.resize(n); vl r(n); if (n) r[0] = 1;
	for (; b; b /= 2, a = nMulT(a, a, n)) if (b & 1) r = nMulT(r, a, n);
	return r;
}
vl nRem(vl f, const vl& g) { // g.back() != 0
	ll il = minv(g.back()); int m = sz(g);
	for (int i = sz(f) - 1; i >= m - 1; i--) {
		ll q = f[i] * il % M;
		rep(j,0,m) f[i - m + 1 + j] = ((f[i - m + 1 + j] - q * g[j]) % M + M) % M;
	}
	if (sz(f) > m - 1) f.resize(m - 1);
	return f;
}
vl nQuo(vl f, const vl& g) {
	ll il = minv(g.back()); int m = sz(g);
	vl q(max(0, sz(f) - m + 1));
	for (int i = sz(f) - 1; i >= m - 1; i--) {
		ll c = q[i - m + 1] = f[i] * il % M;
		rep(j,0,m) f[i - m + 1 + j] = ((f[i - m + 1 + j] - c * g[j]) % M + M) % M;
	}
	return q;
}
bool eqTrim(vl a, vl b) {
	while (sz(a) && !a.back()) a.pop_back();
	while (sz(b) && !b.back()) b.pop_back();
	return a == b;
}

void testBasic() {
	rep(it,0,20000) {
		int n = (int)rnd(0, 8), m = (int)rnd(0, 8);
		vl a = randPoly(n), b = randPoly(m);
		poly pa = P(a), pb_ = P(b);
		ll c = coef(), x = coef();
		assert(V(rev(pa)) == vl(a.rbegin(), a.rend()));
		{ vl t = a; while (sz(t) && !t.back()) t.pop_back();
		  assert(V(REMZ(pa)) == t); poly q = pa; remz(q); assert(V(q) == t); }
		{ int s = (int)rnd(0, 5); vl t(s, 0); t.insert(t.end(), all(a)); assert(V(shift(pa, s)) == t);
		  s = (int)rnd(0, n); assert(V(shift(pa, -s)) == vl(a.begin() + s, a.end())); }
		{ int s = (int)rnd(0, 12); vl t = a; t.resize(s); assert(V(RSZ(pa, s)) == t); }
		assert(eval(pa, mint(x)).v == nEval(a, x));
		{ vl d; rep(i,1,n) d.pb(i * a[i] % M); assert(V(dif(pa)) == d);
		  vl g(n + 1); rep(i,0,n) g[i + 1] = a[i] * minv(i + 1) % M; assert(V(integ(pa)) == g);
		  assert(V(dif(integ(pa))) == a); }
		vl s(max(n, m)), d(max(n, m));
		rep(i,0,max(n, m)) {
			ll u = i < n ? a[i] : 0, v = i < m ? b[i] : 0;
			s[i] = (u + v) % M, d[i] = (u - v + M) % M;
		}
		assert(V(pa + pb_) == s && V(pa - pb_) == d);
		{ poly q = pa; q += pb_; assert(V(q) == s); q = pa; q -= pb_; assert(V(q) == d); }
		{ vl t = a; for (ll& y : t) y = (M - y) % M; assert(V(-pa) == t); }
		{ vl t = a; for (ll& y : t) y = y * c % M;
		  assert(V(pa * mint(c)) == t && V(mint(c) * pa) == t);
		  poly q = pa; q *= mint(c); assert(V(q) == t);
		  if (c) { assert(V(q / mint(c)) == a); q /= mint(c); assert(V(q) == a); } }
		assert(V(pa * pb_) == nMul(a, b));
		{ poly q = pa; q *= pb_; assert(V(q) == nMul(a, b)); }
		assert(V(conv(pa, pb_)) == nMul(a, b));
	}
	// integ's inverse table must also work after being grown several times
	for (int n : {1, 50, 3, 2000, 100}) {
		vl a = randPoly(n), g(n + 1); rep(i,0,n) g[i + 1] = a[i] * minv(i + 1) % M;
		assert(V(integ(P(a))) == g);
	}
}

void testFft() {
	// fft matches the DFT definition with w = RT^((MOD-1)/n), and inverts
	assert(mpow(RT, (M - 1) / 2) == M - 1); // RT generates the 2-power part
	for (int n = 1; n <= 64; n *= 2) rep(it,0,30) {
		vl a = randPoly(n); poly p = P(a);
		fft(p);
		ll w = mpow(RT, (M - 1) / n);
		rep(k,0,n) { ll s = 0; rep(j,0,n) s = (s + a[j] * mpow(w, (ll)j * k)) % M; assert(p[k].v == s); }
		fft(p, 1); assert(V(p) == a);
	}
	rep(it,0,300) {
		vl a = randPoly((int)rnd(1, 300)), b = randPoly((int)rnd(1, 300));
		assert(V(conv(P(a), P(b))) == nMul(a, b));
	}
	{ // all coefficients MOD-1, power-of-two-boundary sizes
		for (int n : {1, 2, 3, 4, 5, 127, 128, 129, 1024, 1025}) {
			vl a(n, M - 1), b(n + 1, M - 1);
			assert(V(conv(P(a), P(b))) == nMul(a, b));
		}
	}
}

void testSeries() {
	rep(it,0,6000) {
		int n = (int)rnd(1, it < 5000 ? 20 : 150), m = (int)rnd(1, 2 * n);
		vl a = randPoly(m);
		// inv
		a[0] = rnd(1, M - 1);
		assert(V(inv(P(a), n)) == nInv(a, n));
		{ vl t = nMulT(a, V(inv(P(a), n)), n), one(n); one[0] = 1; assert(t == one); }
		// sqrt / log need A[0] = 1
		a[0] = 1;
		{ vl b = V(sqrt(P(a), n)); assert(sz(b) == n && b[0] == 1);
		  vl t = a; t.resize(n); assert(nMulT(b, b, n) == t); }
		assert(V(log(P(a), n)) == nLog(a, n));
		// exp needs A[0] = 0
		a[0] = 0;
		assert(V(exp(P(a), n)) == nExp(a, n));
		{ vl e = nExp(a, n); e.resize(n); vl l = nLog(e, n), aa = a; aa.resize(n); assert(l == aa); }
	}
}

void testPow() {
	rep(it,0,6000) {
		int n = (int)rnd(1, 24);
		int m = (int)rnd(1, 2 * n); // A may be shorter than n
		vl a = randPoly(m);
		int z = (int)rnd(0, 3) ? 0 : (int)rnd(0, m);
		rep(i,0,z) a[i] = 0;
		ll b;
		switch (rnd(0, 5)) {
			case 0: b = rnd(0, 3); break;
			case 1: b = rnd(0, 2 * n); break;
			case 2: b = M + rnd(-2, 2); break;
			case 3: b = (M - 1) * rnd(1, 3) + rnd(-1, 1); break;
			case 4: b = rnd(0, (ll)1e18); break;
			default: b = M * rnd(1, 1000) + rnd(0, 3);
		}
		vl got = V(pow(P(a), b, n)), want = nPow(a, b, n);
		if (got != want) {
			cerr << "pow mismatch n=" << n << " b=" << b << " a="; for (ll x : a) cerr << x << ' ';
			cerr << "\n got "; for (ll x : got) cerr << x << ' ';
			cerr << "\nwant "; for (ll x : want) cerr << x << ' ';
			cerr << endl; assert(0);
		}
	}
}

void testDivision() {
	rep(it,0,20000) {
		int m = (int)rnd(1, 12), n = (int)rnd(0, 30);
		vl f = randPoly(n), g = randPoly(m);
		g.back() = rnd(1, M - 1); // divisor needs a non-zero leading coefficient
		auto qr = quoRem(P(f), P(g));
		vl q = V(qr.first), r = V(qr.second);
		if (n < m) assert(q.empty() && r == f);
		else {
			assert(q == nQuo(f, g) && r == nRem(f, g));
			assert(sz(q) == n - m + 1 && sz(r) == m - 1);
		}
		assert(V(mod(P(f), P(g))) == r);
	}
	// x^k mod f
	rep(it,0,3000) {
		int m = (int)rnd(2, 10);
		vl f = randPoly(m); f.back() = rnd(1, M - 1);
		ll k = rnd(0, 3) ? rnd(0, 40) : rnd(0, (ll)1e18);
		vl want;
		if (k <= 40) { vl xk(k + 1); xk[k] = 1; want = nRem(xk, f); }
		else {
			vl r{1}, a{0, 1};
			for (ll e = k; e; e /= 2, a = nRem(nMul(a, a), f)) if (e & 1) r = nRem(nMul(r, a), f);
			want = r;
		}
		assert(eqTrim(V(xkmodf(k, P(f))), want));
	}
#ifdef QUOREM_LEADING_ZERO
	{ // divisor with a zero leading coefficient: (x^2 + 1) / (x + 0x^2)
		auto qr = quoRem(P({1, 0, 1}), P({0, 1, 0}));
		assert(eqTrim(V(qr.first), {0, 1}) && eqTrim(V(qr.second), {1}));
	}
#endif
	// linear recurrences, c is 1-indexed: a[k] = c[1]a[k-1] + ... + c[n]a[k-n]
	rep(it,0,3000) {
		int n = (int)rnd(1, 8);
		vl s = randPoly(n), c = randPoly(n + 1);
		ll k = rnd(0, 60);
		vl a = s;
		rep(i,n,(int)k+1) { ll t = 0; rep(j,1,n+1) t = (t + c[j] * a[i - j]) % M; a.pb(t); }
		assert(solve_linrec(P(s), P(c), n, k).v == a[k]);
	}
	{ // Fibonacci with a large index, checked by 2x2 matrix doubling
		ll k = (ll)1e18, a = 0, b = 1; // (F(i), F(i+1))
		for (int i = 62; i >= 0; i--) {
			ll c = a * ((2 * b - a + M) % M) % M, d = (a * a + b * b) % M;
			a = c, b = d;
			if (k >> i & 1) { ll t = (a + b) % M; a = b, b = t; }
		}
		assert(solve_linrec(P({0, 1}), P({0, 1, 1}), 2, k).v == a);
	}
}

void testMultipoint() {
	rep(it,0,4000) {
		int n = (int)rnd(0, 25), m = (int)rnd(1, 25);
		vl f = randPoly(n), x(m);
		for (ll& y : x) y = rnd(0, 3) ? rnd(0, M - 1) : rnd(0, 3); // duplicates allowed
		vl got = V(multiEval(P(f), P(x)));
		assert(sz(got) == m);
		rep(i,0,m) assert(got[i] == nEval(f, x[i]));
	}
	assert(multiEval(P({1, 2, 3}), {}).empty());
	assert(interpolate(vector<pair<T, T>>{}).empty());
	assert(V(interpolate({{mint(5), mint(7)}})) == vl{7});
	rep(it,0,4000) {
		int n = (int)rnd(1, 25);
		set<ll> xs;
		while (sz(xs) < n) xs.insert(rnd(0, 1) ? rnd(0, M - 1) : rnd(0, 40));
		vl x(all(xs)); shuffle(all(x), rng);
		vl f = randPoly(n); // unique polynomial of degree < n through the points
		vector<pair<T, T>> pts;
		rep(i,0,n) pts.pb({mint(x[i]), mint(nEval(f, x[i]))});
		vl got = V(interpolate(pts));
		assert(sz(got) <= n && eqTrim(got, f));
	}
	for (int n : {1000, 4096, 5000}) {
		vl f = randPoly(n), x(n);
		iota(all(x), rnd(0, M - 1 - n));
		vl got = V(multiEval(P(f), P(x)));
		vector<pair<T, T>> pts;
		rep(i,0,n) { assert(got[i] == nEval(f, x[i])); pts.pb({mint(x[i]), mint(got[i])}); }
		assert(eqTrim(V(interpolate(pts)), f));
	}
}

void testLarge() {
	// identities at sizes where the naive oracle is too slow
	int n = 1 << 16;
	vl a = randPoly(n + 5); a[0] = 1;
	poly A = P(a), one(n); one[0] = 1;
	assert(V(RSZ(conv(A, inv(A, n)), n)) == V(one));
	poly s = sqrt(A, n);
	assert(V(RSZ(conv(s, s), n)) == V(RSZ(A, n)));
	assert(V(exp(log(A, n), n)) == V(RSZ(A, n)));
	poly e = pow(A, 5, n), a2 = RSZ(conv(A, A), n), a4 = RSZ(conv(a2, a2), n);
	assert(V(e) == V(RSZ(conv(a4, A), n)));
	{ // conv against a random-point evaluation check
		poly B = P(randPoly(n)), C = conv(A, B);
		rep(it,0,5) { mint x(rnd(0, M - 1)); assert((eval(A, x) * eval(B, x)).v == eval(C, x).v); }
		auto qr = quoRem(C, B); // exact division
		assert(V(qr.first) == V(A) && V(qr.second) == vl(n - 1, 0));
	}
}

int main() {
	testBasic();
	testFft();
	testSeries();
	testPow();
	testDivision();
	testMultipoint();
	testLarge();
	cout << "Tests passed!" << endl;
}
