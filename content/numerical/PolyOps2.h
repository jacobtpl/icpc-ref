/**
 * Author: jacobtpl, Benq
 * Date: 2024
 * Description: Operations on formal power series
 */
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
poly inv(poly A, int n) { // Q-(1/Q-A)/(-Q^{-2})
	poly B{invert(A[0])};
	for (int x = 2; x/2 < n; x *= 2)
		B = 2*B-RSZ(conv(RSZ(A,x),conv(B,B)),x);
	return RSZ(B,n);
} // 455ce0
poly sqrt(const poly& A, int n) {  // Q-(Q^2-A)/(2Q)
	assert(A[0].v == 1); poly B{1};
	for (int x = 2; x/2 < n; x *= 2)
		B = invert(T(2))*RSZ(B+conv(RSZ(A,x),inv(B,x)),x);
	return RSZ(B,n);
} // ecb486
// return {quotient, remainder}
pair<poly,poly> quoRem(const poly& f, const poly& g) {
	if (sz(f) < sz(g)) return {{},f};
	poly q = conv(inv(rev(g),sz(f)-sz(g)+1),rev(f));
	q = rev(RSZ(q,sz(f)-sz(g)+1));
	poly r = RSZ(f-conv(q,g),sz(g)-1); return {q,r};
} // bffc5b
poly log(poly A, int n) { assert(A[0].v == 1); // (ln A)' = A'/A
	A.resize(n); return integ(RSZ(conv(dif(A),inv(A,n-1)),n-1));
} // bda418
poly exp(poly A, int n) { assert(A[0].v == 0);
	poly B{1}, IB{1}; // inverse of B
	for (int x = 1; x < n; x *= 2) {
		IB = 2*IB-RSZ(conv(B,conv(IB,IB)),x);
		poly Q = dif(RSZ(A,x)); Q += RSZ(conv(IB,dif(B)-conv(B,Q)),2*x-1); 
		/// first x-1 terms of dif(B)-conv(B,Q) are zero
		B = B+RSZ(conv(B,RSZ(A,2*x)-integ(Q)),2*x); 
	} /// We know that Q=A' is B'/B to x-1 places, we want to find B'/B to 2x-1 places
	return RSZ(B,n);
} // 203953
poly pow(poly A, ll b, int n) {
    if (b==0) { poly r(n,0); r[0]=1; return r; }
    int t = -1;
    for (int i = 0; i < n; i++) if (A[i].v != 0) { t = i; break; }
    if (t == -1) return poly(n, 0);
    mint fac = A[t];
    for (int i = 0; i < n; i++) A[i] /= fac;
    poly p(A.begin()+t, A.end());
    p.resize(n);
    poly q = log(p, n);
    poly r = exp(q * mint(b), n) * pow(fac, b);
    if (t == 0) return r;
    if (b >= n || b*t >= n) return poly(n, 0);
    r.insert(r.begin(), t*b, mint(0));
    r.resize(n);
    return r;
} // 8dc32d
poly mod(const poly& f, const poly& g) { return quoRem(f,g).second; }
poly xkmodf(ll k, poly f) {
    poly r{1}, a{0,1};
    for(;k;k>>=1) {
        if(k&1) r = mod(conv(r,a), f);
        a = mod(conv(a,a), f);
    }
    return r;
} // ef2278
// solve recurrence with initial vals s[0], s[1]... s[n-1]
// a[k] = c[1]*a[k-1] + c[2]*a[k-2] + ... c[n]*a[k-n]
mint solve_linrec(vector<mint> s, vector<mint> c, int n, ll k) {
    poly f(n+1, 0);
    f[n] = 1;
    for (int i=0;i<n;i++) f[i] = mint(-c[n-i]);
    poly r = xkmodf(k, f); r.resize(n);
    mint ans(0);
    for (int i = 0; i < n; i++) ans += r[i] * mint(s[i]);
    return ans;
} // 902e38
