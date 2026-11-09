#include "../utilities/template.h"

#include "../../content/data-structures/SegmentTreeBeats.h"

static unsigned RA = 1231231;
int ra() {
	RA *= 574841;
	RA += 14;
	return (int)(RA >> 1);
}
ll rl() { return ((ll)ra() << 31) ^ ra(); }

const ll INF = LLONG_MAX;

void check(Beats& t, vector<ll>& v, int L, int R) {
	ll ma = -INF, mi = INF, su = 0;
	rep(i,L,R) ma = max(ma, v[i]), mi = min(mi, v[i]), su += v[i];
	assert(t.query(L, R, 0) == ma);
	assert(t.query(L, R, 1) == mi);
	assert(t.query(L, R, 2) == su);
}

// n: array size, A: |a_i| bound, B: |v| bound, q: #ops
void run(int n, ll A, ll B, int q) {
	vector<ll> v(n);
	for (ll& x : v) x = rl() % (2*A + 1) - A;
	Beats t(v);
	rep(it,0,q) {
		int L = ra() % (n + 1), R = ra() % (n + 1);
		if (ra() % 4) { if (L > R) swap(L, R); }
		int k = ra() % 3;
		if (ra() % 3) {
			ll x = rl() % (2*B + 1) - B;
			if (k == 2 && ra() % 2) x = rl() % 3 - 1;
			t.update(L, R, x, k);
			rep(i,L,R) v[i] = k == 0 ? min(v[i], x) :
			                 k == 1 ? max(v[i], x) : v[i] + x;
		} else check(t, v, L, R);
	}
	check(t, v, 0, n);
	rep(i,0,n) check(t, v, i, i + 1);
	assert(t.query(0, 0, 0) == -INF);
	assert(t.query(0, 0, 1) == INF);
	assert(t.query(0, 0, 2) == 0);
}

int main() {
	rep(it,0,20) run(0, 10, 10, 50);
	rep(it,0,500) run(1, 10, 10, 50);
	rep(it,0,3000) run(ra() % 8 + 1, 3, 3, 60);
	rep(it,0,2000) run(ra() % 20 + 1, 2, 2, 100);
	rep(it,0,1000) run(ra() % 20 + 1, 1000000000, 1000000000, 100);
	rep(it,0,300) run(ra() % 100 + 1, 1000000000, 1000000000, 300);
	// large values: |a_i| and |v| near 2^60, sums stay below 2^63
	rep(it,0,500) run(ra() % 6 + 1, 1LL << 60, 1LL << 59, 60);
	rep(it,0,50) run(ra() % 2000 + 1, 1000000000, 1000000000, 1000);
	rep(it,0,3) run(10000, 1000000000000000LL, 1000000000000LL, 5000);
	cout << "Tests passed!" << endl;
}
