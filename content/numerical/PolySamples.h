/**
 * Author: Devin
 * Date: 2026-11-21
 * Description: A polynomial f of degree < n given by its values
 *   y[i] = f(i), i = 0..n-1 (n < MOD, MOD prime).
 *   lagrangeAt needs only ModInt.h; shiftSamples needs PolyOps.h.
 * Time: O(n) for lagrangeAt, O((n+m) log(n+m)) for shiftSamples
 * Status: stress-tested
 */
// w[i] = (-1)^(n-1-i) / (i! (n-1-i)!)
vector<mint> lagrangeWeights(int n) {
	vector<mint> r(n), w(n); mint f = 1; rep(i,1,n) f *= i;
	if (n) r[n-1] = invert(f);
	for (int i = n-1; i > 0; i--) r[i-1] = r[i]*i;
	rep(i,0,n) w[i] = r[i]*r[n-1-i]*((n-1-i)%2 ? -1 : 1);
	return w;
}
// f(x)
mint lagrangeAt(const vector<mint>& y, mint x) {
	int n = sz(y); if (x.v < n) return y[x.v];
	vector<mint> w = lagrangeWeights(n), suf(n+1,1);
	for (int i = n; i--;) suf[i] = suf[i+1]*(x-i);
	mint res = 0, pre = 1;
	rep(i,0,n) res += y[i]*w[i]*pre*suf[i+1], pre *= x-i;
	return res;
}
// f(c), f(c+1), ..., f(c+m-1)
vector<mint> shiftSamples(const vector<mint>& y, mint c, int m) {
	int n = sz(y), k = n+m-1; if (!n || !m) return vector<mint>(m);
	vector<mint> d = lagrangeWeights(n), u(k), pre(k+1,1), ipre(k+1), h(k);
	rep(j,0,k) { // u[j] = c-(n-1)+j, with 0 replaced by 1
		u[j] = c-(n-1)+j; if (!u[j].v) u[j] = 1;
		pre[j+1] = pre[j]*u[j];
	}
	ipre[k] = invert(pre[k]); // ipre[j] = 1/pre[j], h[j] = 1/u[j]
	for (int j = k; j--;) ipre[j] = ipre[j+1]*u[j], h[j] = ipre[j+1]*pre[j];
	rep(i,0,n) d[i] *= y[i];
	d = conv(d,h); vector<mint> res(m);
	rep(j,0,m) { // f(x) = sum_i y[i]w[i]/(x-i) * prod_i (x-i)
		mint x = c+j;
		res[j] = x.v < n ? y[x.v] : d[j+n-1]*pre[j+n]*ipre[j];
	}
	return res;
}
