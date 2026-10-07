#include "../utilities/template.h"

#include "../../content/combinatorial/IntPerm.h"

// independent O(n^2) mixed-radix code: digit i = #{j < i : v[j] > v[i]}
int naive(const vi& v) {
	ll r = 0;
	rep(i,0,sz(v)) {
		int c = 0;
		rep(j,0,i) c += v[j] > v[i];
		r = r * (i + 1) + c;
	}
	assert(r <= INT_MAX);
	return (int)r;
}

int main() {
	mt19937 rng(12345);
	vi e;
	assert(permToInt(e) == 0);
	// exhaustive: bijection onto [0, n!)
	int fact = 1;
	rep(n,1,11) {
		fact *= n;
		vi v(n); iota(all(v), 0);
		vector<bool> seen(fact);
		int cnt = 0;
		do {
			int r = permToInt(v);
			assert(0 <= r && r < fact);
			assert(!seen[r]); seen[r] = 1;
			if (n <= 8 || cnt % 97 == 0) assert(r == naive(v));
			cnt++;
		} while (next_permutation(all(v)));
		assert(cnt == fact);
	}
	// random, n = 11, 12 (12! - 1 is the largest value that fits in int)
	rep(n,11,13) {
		ll f = 1; rep(i,1,n+1) f *= i;
		vi v(n); iota(all(v), 0);
		assert(permToInt(v) == 0);
		reverse(all(v));
		assert(permToInt(v) == f - 1);
		set<pair<int, vi>> seen;
		rep(it,0,200000) {
			shuffle(all(v), rng);
			int r = permToInt(v);
			assert(0 <= r && r < f && r == naive(v));
		}
	}
	cout << "Tests passed!" << endl;
}
