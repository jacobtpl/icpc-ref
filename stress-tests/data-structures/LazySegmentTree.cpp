#include "../utilities/template.h"

#include "../../content/data-structures/LazySegmentTree.h"

static unsigned R;
int ra() {
	R *= 791231;
	R += 1231;
	return (int)(R >> 1);
}

volatile int res;
int main() {
	int N = 10;
	vi v(N);
	iota(all(v), 0);
	random_shuffle(all(v), [](int x) { return ra() % x; });
	Node* tr = new Node(v,0,N);
	rep(i,0,N) rep(j,0,N) if (i <= j) {
		int ma = -inf;
		rep(k,i,j) ma = max(ma, v[k]);
		assert(ma == tr->query(i,j));
	}
	rep(it,0,1000000) {
		int i = ra() % (N+1), j = ra() % (N+1);
		if (i > j) swap(i, j);
		int x = (ra() % 10) - 5;

		int r = ra() % 100;
		if (r < 30) {
			::res = tr->query(i, j);
			int ma = -inf;
			rep(k,i,j) ma = max(ma, v[k]);
			assert(ma == ::res);
		}
		else if (r < 70) {
			tr->add(i, j, x);
			rep(k,i,j) v[k] += x;
		}
		else {
			tr->set(i, j, x);
			rep(k,i,j) v[k] = x;
		}
	}
	// random sizes and offsets, both constructors, larger values
	mt19937 rng(3);
	rep(it,0,6000) {
		int n = (int)(rng() % 40) + 1, lo = 0;
		if (it % 100 == 0) n = 2000;
		bool sparse = it % 2;
		int V = it % 3 == 0 ? 5 : 1000000;
		vi w(n, -inf);
		Node* t;
		if (sparse) {
			if (it % 4 == 1) lo = (int)(rng() % 2000001) - 1000000;
			else if (it % 8 == 3) lo = INT_MAX - n; // near the int limit
			else if (it % 8 == 7) lo = INT_MIN;
			t = new Node(lo, lo + n);
		} else {
			for (int& x : w) x = (int)(rng() % (2 * V + 1)) - V;
			t = new Node(w, 0, n);
		}
		rep(q,0,n > 100 ? 4000 : 120) {
			int i = (int)(rng() % (n + 1)), j = (int)(rng() % (n + 1));
			if (i > j) swap(i, j);
			if (rng() % 10 == 0) i = 0, j = n;
			int x = (int)(rng() % (2 * V + 1)) - V, r = (int)(rng() % 3);
			if (r == 0) {
				int ma = -inf;
				rep(k,i,j) ma = max(ma, w[k]);
				assert(t->query(lo + i, lo + j) == ma);
			} else if (r == 1) {
				x = x % 1000; // keep sums far from overflow
				// untouched sparse cells are -inf; values must never go below that
				if (i < j && *min_element(w.begin() + i, w.begin() + j) + x < -inf) x = -x;
				t->add(lo + i, lo + j, x);
				rep(k,i,j) w[k] += x;
			} else {
				t->set(lo + i, lo + j, x);
				rep(k,i,j) w[k] = x;
			}
		}
		rep(i,0,n) assert(t->query(lo + i, lo + i + 1) == w[i]);
		// ranges sticking out of [lo, hi) are clipped
		if (!sparse) assert(t->query(-5, n + 5) == *max_element(all(w)));
	}
	cout<<"Tests passed!"<<endl;
}
