/**
 * Author: jacobtpl, Benq
 * Date: 2024
 * Description: Operations on formal power series mod 998244353.
 *   Needs ModInt.h. A poly is its coefficient list, lowest degree first;
 *   f(A, n) returns the first n coefficients of f(A).
 *   Whole file is O(n log n) per operation unless noted.
 */
int const RT = 5;
using T = mint;
using poly = vector<mint>;
void remz(poly& p) { while (sz(p)&&p.back().v==0) p.pop_back(); }
poly REMZ(poly p) { remz(p); return p; }
poly rev(poly p) { reverse(all(p)); return p; }
poly shift(poly p, int x) { 
	if (x >= 0) p.insert(begin(p),x,0); 
	else assert(sz(p)+x >= 0), p.erase(begin(p),begin(p)-x);
	return p; 
}
poly RSZ(const poly& p, int x) { 
	if (x <= sz(p)) return poly(begin(p),begin(p)+x);
	poly q = p; q.resize(x); return q;  }
T eval(const poly& p, T x) { // evaluate at point x
	T res = 0; for (int i = sz(p)-1; i>=0; i--) res = x*res+p[i]; 
	return res; }
poly dif(const poly& p) { // differentiate
	poly res; rep(i,1,sz(p)) res.pb(T(i)*p[i]); 
	return res; }
poly integ(const poly& p) { // integrate
	static poly invs{0,1};
	for (int i = sz(invs); i <= sz(p); ++i) 
		invs.pb(-MOD/i*invs[MOD%i]);
	poly res(sz(p)+1); rep(i,0,sz(p)) res[i+1] = p[i]*invs[i+1];
	return res; 
}
poly& operator+=(poly& l, const poly& r) {
	l.resize(max(sz(l),sz(r))); rep(i,0,sz(r)) l[i] += r[i]; 
	return l; }
poly& operator-=(poly& l, const poly& r) {
	l.resize(max(sz(l),sz(r))); rep(i,0,sz(r)) l[i] -= r[i]; 
	return l; }
poly& operator*=(poly& l, const T& r) { for (auto &t:l) t *= r; 
	return l;	 }
poly& operator/=(poly& l, const T& r) { for (auto &t:l) t /= r; 
	return l; }
poly operator+(poly l, const poly& r) { return l += r; }
poly operator-(poly l, const poly& r) { return l -= r; }
poly operator-(poly l) { for (auto &t:l) t *= -1; return l; }
poly operator*(poly l, const T& r) { return l *= r; }
poly operator*(const T& r, const poly& l) { return l*r; }
poly operator/(poly l, const T& r) { return l /= r;	 } // 8d2d0d
void fft(vector<T>& A, bool inverse = 0) { // NTT
	int n = sz(A); assert((MOD-1)%n == 0); vector<T> B(n);
	for(int b = n/2; b; b /= 2, swap(A,B)) { // w = n/b'th root
		T w = pow(mint(RT),(MOD-1)/n*b), m = 1; 
		for(int i = 0; i < n; i += b*2, m *= w) rep(j,0,b) {
			T u = A[i+j], v = A[i+j+b]*m;
			B[i/2+j] = u+v; B[i/2+j+n/2] = u-v;
		}
	}
	if (inverse) { reverse(1+all(A)); 
		T z = invert(T(n)); for (auto &t:A) t *= z; }
} // for NTT-able moduli -- 3397f7
vector<T> conv(vector<T> A, vector<T> B) {
	if (!min(sz(A),sz(B))) return {};
	int s = sz(A)+sz(B)-1, n = 1; for (; n < s; n *= 2);
	A.resize(n), fft(A); B.resize(n), fft(B);
	rep(i,0,n) A[i] *= B[i];
	fft(A,1); A.resize(s); return A;
} // 7956e1
poly operator*(const poly& l, const poly& r) {
	assert(MOD == 998244353); return conv(l,r);
} // b6ba5b
// poly operator*(const poly& l, const poly& r) { // O(NM)
// 	if (!min(sz(l),sz(r))) return {};
// 	poly x(sz(l)+sz(r)-1);
// 	rep(i,0,sz(l)) rep(j,0,sz(r)) x[i+j] += l[i]*r[j];
// 	return x;
// }
// Other prime MOD (e.g. 1e9+7): replace conv by this, which needs
// FastFourierTransform.h and FastFourierTransformMod.h, and drop the
// assert in operator*. Everything below then works unchanged.
// vector<T> conv(const vector<T>& A, const vector<T>& B) {
// 	vl a, b; for (T x : A) a.pb(x.v); for (T x : B) b.pb(x.v);
// 	vl c = convMod<MOD>(a,b); return poly(all(c));
// }
poly& operator*=(poly& l, const poly& r) { return l = l*r; }
poly inv(poly A, int n) { // needs A[0] != 0; Q-(1/Q-A)/(-Q^{-2})
	poly B{invert(A[0])};
	for (int x = 2; x/2 < n; x *= 2)
		B = 2*B-RSZ(conv(RSZ(A,x),conv(B,B)),x);
	return RSZ(B,n);
} // 455ce0
T sqrt(T a) { // a root of a != 0, or 0 if there is none
	if (pow(a,(MOD-1)/2).v != 1) return 0;
	int s = __builtin_ctz(MOD-1), m; T z = 2;
	while (pow(z,(MOD-1)/2).v == 1) z += 1;
	T c = pow(z,MOD>>s), x = pow(a,((MOD>>s)+1)/2), t = pow(a,MOD>>s);
	for (; t.v != 1; s = m) {
		T u = t; for (m = 0; u.v != 1; m++) u *= u;
		T b = pow(c,1<<(s-m-1)); x *= b, c = b*b, t *= c;
	}
	return x;
} // 6ab367
poly sqrt(poly A, int n) { // {} if there is none; Q-(Q^2-A)/(2Q)
	int t = 0; while (t < sz(A) && !A[t].v) t++;
	if (t == sz(A)) return poly(n);
	T c = sqrt(A[t]); if (t%2 || !c.v) return {};
	if (t/2 >= n) return poly(n);
	A.erase(begin(A),begin(A)+t); poly B{c};
	for (int x = 2; x/2 < n-t/2; x *= 2)
		B = invert(T(2))*RSZ(B+conv(RSZ(A,x),inv(B,x)),x);
	B.resize(n-t/2); return shift(B,t/2);
} // a59092
// return {quotient, remainder}; needs g.back() != 0
pair<poly,poly> quoRem(const poly& f, const poly& g) {
	if (sz(f) < sz(g)) return {{},f};
	poly q = conv(inv(rev(g),sz(f)-sz(g)+1),rev(f));
	q = rev(RSZ(q,sz(f)-sz(g)+1));
	poly r = RSZ(f-conv(q,g),sz(g)-1); return {q,r};
} // bffc5b
poly log(poly A, int n) { // needs A[0] == 1; (ln A)' = A'/A
	if (!n) return {};
	assert(A[0].v == 1); A.resize(n);
	return integ(RSZ(conv(dif(A),inv(A,n-1)),n-1));
} // 970f4a
poly exp(poly A, int n) { assert(A.empty() || A[0].v == 0);
	poly B{1}, IB{1}; // inverse of B
	for (int x = 1; x < n; x *= 2) {
		IB = 2*IB-RSZ(conv(B,conv(IB,IB)),x);
		poly Q = dif(RSZ(A,x)); Q += RSZ(conv(IB,dif(B)-conv(B,Q)),2*x-1); 
		/// first x-1 terms of dif(B)-conv(B,Q) are zero
		B = B+RSZ(conv(B,RSZ(A,2*x)-integ(Q)),2*x); 
	} /// We know that Q=A' is B'/B to x-1 places, we want to find B'/B to 2x-1 places
	return RSZ(B,n);
} // 07e710
poly pow(poly A, ll b, int n) { // A^b mod x^n, any A, b >= 0
	if (!n) return {};
	A.resize(n); int t = 0; while (t < n && !A[t].v) t++;
	if (!b) { poly r(n); r[0] = 1; return r; }
	if (t == n || (t && b >= n) || t*b >= n) return poly(n);
	int k = int(n-t*b); mint c = A[t];
	poly p = RSZ(poly(begin(A)+t, end(A))/c, k);
	p = exp(log(p,k)*mint(b),k)*pow(c,b);
	p.insert(begin(p), t*b, 0); return p;
} // 0f5d41
poly mod(const poly& f, const poly& g) { return quoRem(f,g).second; }
poly xkmodf(ll k, poly f) { // x^k mod f
    poly r{1}, a{0,1};
    for(;k;k>>=1) {
        if(k&1) r = mod(conv(r,a), f);
        a = mod(conv(a,a), f);
    }
    return r;
} // ef2278
// [x^k] P/Q in O(d log d log k), Q[0] != 0, deg P < deg Q = d
T kthCoef(poly P, poly Q, ll k) {
	for (; k; k /= 2) {
		poly R = Q; for (int i = 1; i < sz(R); i += 2) R[i] = -R[i];
		poly A = conv(P,R), B = conv(Q,R); P.clear(), Q.clear();
		for (int i = k&1; i < sz(A); i += 2) P.pb(A[i]);
		for (int i = 0; i < sz(B); i += 2) Q.pb(B[i]);
	}
	return sz(P) ? P[0]/Q[0] : 0;
} // 993413
// solve recurrence with initial vals s[0], s[1]... s[n-1]
// a[k] = c[1]*a[k-1] + c[2]*a[k-2] + ... c[n]*a[k-n]
mint solve_linrec(vector<mint> s, vector<mint> c, int n, ll k) {
	poly Q(n+1); Q[0] = 1; rep(i,1,n+1) Q[i] = -c[i];
	return kthCoef(RSZ(conv(RSZ(s,n),Q),n),Q,k);
} // dc8033
// shortest c (c[0] = 0) with s[i] = sum of c[j]*s[i-j] over j >= 1
// for all i >= sz(c)-1, O(n^2). Feed at least 2x the order, then
// solve_linrec(s, c, sz(c)-1, k) is s[k].
poly berlekampMassey(const poly& s) {
	int n = sz(s), L = 0, m = 0; poly C(n+1), B(n+1), t;
	C[0] = B[0] = 1; T b = 1;
	rep(i,0,n) { ++m;
		T d = s[i]; rep(j,1,L+1) d += C[j]*s[i-j];
		if (!d.v) continue;
		t = C; T q = d/b; rep(j,m,n+1) C[j] -= q*B[j-m];
		if (2*L > i) continue;
		L = i+1-L, B = t, b = d, m = 0;
	}
	C.resize(L+1); C[0] = 0; return -C;
} // c1cead
// returns f(x+c)
poly taylorShift(poly f, T c) {
	int n = sz(f); poly F(n+1,1), A(n), B(n); T p = 1;
	rep(i,1,n+1) F[i] = F[i-1]*T(i);
	rep(i,0,n) A[n-1-i] = f[i]*F[i], B[i] = p/F[i], p *= c;
	A = conv(A,B); rep(i,0,n) f[i] = A[n-1-i]/F[i];
	return f;
} // 293ba2
