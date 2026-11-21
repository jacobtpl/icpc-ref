#include "../utilities/template.h"

#include "../../content/various/MaxPlusConvolution.h"

vector<ll> brute(const vector<ll>& a, const vector<ll>& b) {
	vector<ll> c(sz(a) + sz(b) - 1, LLONG_MIN);
	rep(i,0,sz(a)) rep(j,0,sz(b)) c[i+j] = max(c[i+j], a[i] + b[j]);
	return c;
}
vector<ll> concave(int n, int w) {
	vector<ll> d(n);
	for (auto &x : d) x = rand() % (2 * w + 1) - w;
	sort(d.rbegin(), d.rend());
	d[0] = rand() % 21 - 10;
	rep(i,1,n) d[i] += d[i-1];
	return d;
}
vector<ll> arb(int n, int w) {
	vector<ll> d(n);
	for (auto &x : d) x = rand() % (2 * w + 1) - w;
	return d;
}

int main() {
	rep(it,0,300000) {
		int n = rand() % 12 + 1, m = rand() % 12 + 1;
		int w = it % 3 == 0 ? 1 : it % 3 == 1 ? 3 : 1000;
		auto a = concave(n, w), b = concave(m, w), r = arb(m, w);
		auto c = maxPlusCC(a, b);
		assert(c == brute(a, b));
		rep(i,2,sz(c)) assert(c[i] - c[i-1] <= c[i-1] - c[i-2]);
		assert(maxPlusCA(a, b) == c);
		assert(maxPlusCA(a, r) == brute(a, r));
	}
	rep(it,0,200) {
		int n = rand() % 300 + 1, m = rand() % 300 + 1;
		auto a = concave(n, 5), r = arb(m, 5);
		assert(maxPlusCA(a, r) == brute(a, r));
	}
	cout<<"Tests passed!"<<endl;
}
