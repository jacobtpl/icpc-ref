#include "../utilities/template.h"

#include "../../content/data-structures/FenwickTree.h"

int main() {
	rep(it,0,100000) {
		int N = rand() % 10;
		FT fw(N);
		vi t(N);
		rep(i,0,N) {
			int v = rand() % 3;
			fw.update(i, v);
			t[i] += v;
		}
		int q = rand() % 20;
		int ind = fw.lower_bound(q);
		int res = -1, sum = 0;
		rep(i,0,N+1) {
			if (sum < q) res = i;
			if (i != N) sum += t[i];
		}
		assert(res == ind);
	}
	// update / query against a plain array, with negative and huge values
	mt19937_64 rng(7);
	{ FT e(0); assert(e.query(0) == 0 && e.lower_bound(1) == 0 && e.lower_bound(0) == -1); }
	rep(it,0,30000) {
		int N = (int)(rng() % 12) + (it % 50 == 0 ? 500 : 0);
		ll V = it % 3 == 0 ? 3 : it % 3 == 1 ? 1000000007 : (ll)4e18 / (N + 1) / 40;
		FT fw(N);
		vector<ll> t(N);
		rep(q,0,40) {
			if (N && rng() % 2) {
				int i = (int)(rng() % N);
				ll d = (ll)(rng() % (2 * V + 1)) - V;
				fw.update(i, d), t[i] += d;
			} else {
				int pos = (int)(rng() % (N + 1));
				ll sum = 0;
				rep(i,0,pos) sum += t[i];
				assert(fw.query(pos) == sum);
			}
		}
		ll sum = 0;
		rep(i,0,N) assert(fw.query(i) == sum), sum += t[i];
		assert(fw.query(N) == sum);
	}
	// lower_bound with non-negative values (required), all positions and sums
	rep(it,0,30000) {
		int N = (int)(rng() % 12) + (it % 50 == 0 ? 300 : 0);
		ll V = it % 2 ? 4 : (ll)1e15;
		FT fw(N);
		vector<ll> t(N), pre(N + 1);
		rep(i,0,N) if (rng() % 3) t[i] = (ll)(rng() % V), fw.update(i, t[i]);
		rep(i,0,N) pre[i+1] = pre[i] + t[i];
		auto brute = [&](ll q) {
			if (q <= 0) return -1;
			rep(i,0,N) if (pre[i+1] >= q) return i;
			return N;
		};
		rep(i,0,N+1) for (ll d : {-1, 0, 1}) assert(fw.lower_bound(pre[i] + d) == brute(pre[i] + d));
		for (ll q : {LLONG_MIN, -1LL, 0LL, 1LL, LLONG_MAX}) assert(fw.lower_bound(q) == brute(q));
	}
	{ // large, sizes around powers of two
		for (int N : {1, 2, 3, 1023, 1024, 1025, 1 << 20, (1 << 20) + 1, 1000000}) {
			FT fw(N);
			vector<ll> t(N), pre(N + 1);
			rep(i,0,min(N, 200000)) {
				int j = (int)(rng() % N); ll d = (ll)(rng() % 1000);
				fw.update(j, d), t[j] += d;
			}
			rep(i,0,N) pre[i+1] = pre[i] + t[i];
			rep(i,0,200000) {
				int j = (int)(rng() % (N + 1));
				assert(fw.query(j) == pre[j]);
				ll q = (ll)(rng() % (pre[N] + 5)) - 1;
				int r = fw.lower_bound(q);
				if (q <= 0) assert(r == -1);
				else if (q > pre[N]) assert(r == N);
				else assert(pre[r+1] >= q && pre[r] < q);
			}
		}
	}
	cout<<"Tests passed!"<<endl;
}
