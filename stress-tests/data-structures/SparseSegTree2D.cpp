#include "../utilities/template.h"
#define pb push_back

#include "../../content/data-structures/SparseSegTree2D.h"

mt19937_64 rng(12345);
ll rnd(ll a, ll b) { return a + (ll)(rng() % (unsigned long long)(b - a + 1)); }

// picks coordinates that collide often and hug the borders/midpoints
int coord(int n, const vi& pool) {
	int t = (int)rnd(0, 9);
	if (t == 0) return 0;
	if (t == 1) return n - 1;
	if (t == 2) return (int)rnd(0, n - 1);
	int x = pool[rnd(0, sz(pool) - 1)] + (int)rnd(-1, 1);
	return min(max(x, 0), n - 1);
}

void run(int R, int C, int ops, ll maxv) {
	SegTree2D st(R, C);
	map<pii, ll> m;
	vi pr, pc;
	rep(i,0,4) pr.pb((int)rnd(0, R - 1)), pc.pb((int)rnd(0, C - 1));
	pr.pb(R / 2); pc.pb(C / 2);
	rep(it,0,ops) {
		if (rnd(0, 1)) {
			int p = coord(R, pr), q = coord(C, pc);
			ll k = rnd(0, 3) ? rnd(0, maxv) : 0; // overwrites, also with 0
			st.update(p, q, k);
			m[{p, q}] = k;
		} else {
			int p = coord(R, pr), u = coord(R, pr);
			int q = coord(C, pc), v = coord(C, pc);
			if (rnd(0, 5) == 0) p = 0, u = R - 1;
			if (rnd(0, 5) == 0) q = 0, v = C - 1;
			if (p > u) swap(p, u);
			if (q > v) swap(q, v);
			ll exp = DEFAULT;
			for (auto& [key, val] : m)
				if (p <= key.first && key.first <= u && q <= key.second && key.second <= v)
					exp = func(exp, val);
			ll got = st.query(p, q, u, v);
			if (got != exp) {
				cerr << "R=" << R << " C=" << C << " query(" << p << "," << q << "," << u << "," << v
					<< ") got " << got << " expected " << exp << endl;
				abort();
			}
		}
	}
}

int main() {
	{ // empty structure
		SegTree2D st(1000000000, 1000000000);
		assert(st.query(0, 0, 999999999, 999999999) == 0);
		SegTree2D s1(1, 1);
		assert(s1.query(0, 0, 0, 0) == 0);
		s1.update(0, 0, 7);
		assert(s1.query(0, 0, 0, 0) == 7);
		s1.update(0, 0, 3);
		assert(s1.query(0, 0, 0, 0) == 3);
	}
	const int B = 1000000000;
	vector<pii> dims = {{1,1},{1,2},{2,1},{2,2},{1,7},{7,1},{3,4},{5,5},{8,8},{13,9},
		{B,B},{B,1},{1,B},{1000,B},{B,1000},{B-1,B-3},{1<<30,1<<30}};
	for (auto [R, C] : dims) rep(it,0,300) {
		run(R, C, 60, it % 3 == 0 ? 3 : it % 3 == 1 ? 1000 : (ll)4e18);
	}
	rep(it,0,10) run(B, B, 3000, (ll)4e18);
	rep(it,0,2000) {
		int R = (int)rnd(1, 6), C = (int)rnd(1, 6);
		run(R, C, 80, 5);
	}
	cout<<"Tests passed!"<<endl;
}
