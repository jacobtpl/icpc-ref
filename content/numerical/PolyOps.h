/**
 * Author: jacobtpl, Benq
 * Date: 2024
 * Description: Operations on formal power series
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
poly operator/(poly l, const T& r) { return l /= r;	 }
poly operator*(const poly& l, const poly& r) {
	if (!min(sz(l),sz(r))) return {};
	poly x(sz(l)+sz(r)-1); 
	rep(i,0,sz(l)) rep(j,0,sz(r)) x[i+j] += l[i]*r[j];
	return x;
}
poly& operator*=(poly& l, const poly& r) { return l = l*r; }
