/**
 * Author: Benq
 * Description: Multipoint evaluation and interpolation
 * Time: O(N\log^2 N)
 */
#include "PolyOps.h"

void segProd(vector<poly>& stor, poly& v, int ind, int l, int r) { // v -> places to evaluate at
	if (l == r) { stor[ind] = {mint(0)-v[l],1}; return; }
	int m = (l+r)/2; segProd(stor,v,2*ind,l,m); segProd(stor,v,2*ind+1,m+1,r);
	stor[ind] = conv(stor[2*ind],stor[2*ind+1]);
}
void evalAll(vector<poly>& stor, poly& res, poly v, int ind = 1) {
	v = quoRem(v,stor[ind]).second;
	if (sz(stor[ind]) == 2) { res.pb(sz(v)?v[0]:0); return; }
	evalAll(stor,res,v,2*ind); evalAll(stor,res,v,2*ind+1);
}

// evaluate polynomial v at points in p
poly multiEval(poly v, poly p) {
	vector<poly> stor(4*sz(p)); segProd(stor,p,1,0,sz(p)-1);
	poly res; evalAll(stor,res,v); return res; }

poly combAll(vector<poly>& stor, poly& dems, int ind, int l, int r) {
	if (l == r) return {dems[l]};
	int m = (l+r)/2;
	poly a = combAll(stor,dems,2*ind,l,m), b = combAll(stor,dems,2*ind+1,m+1,r);
	return conv(a,stor[2*ind+1])+conv(b,stor[2*ind]);
}
poly interpolate(vector<pair<T,T>> v) {
	int n = sz(v); poly x; for (auto t:v) x.pb(t.first);
	vector<poly> stor(4*n); segProd(stor,x,1,0,n-1);
	poly dems; evalAll(stor,dems,dif(stor[1]));
	rep(i,0,n) dems[i] = v[i].second/dems[i];
	return combAll(stor,dems,1,0,n-1);
}