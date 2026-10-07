#include "../utilities/template.h"

#include "../../content/data-structures/LineContainer.h"

namespace other {
// source: https://github.com/niklasb/contest-algos/blob/master/convex_hull/dynamic.cpp
const ll is_query = -(1LL<<62);
struct Line {
    ll m, b;
    mutable function<const Line*()> succ;
    bool operator<(const Line& rhs) const {
        if (rhs.b != is_query) return m < rhs.m;
        const Line* s = succ();
        if (!s) return 0;
        ll x = rhs.m;
        return b - s->b < (s->m - m) * x;
    }
};
struct HullDynamic : public multiset<Line> { // will maintain upper hull for maximum
    bool bad(iterator y) {
        auto z = next(y);
        if (y == begin()) {
            if (z == end()) return 0;
            return y->m == z->m && y->b <= z->b;
        }
        auto x = prev(y);
        if (z == end()) return y->m == x->m && y->b <= x->b;
        return (x->b - y->b)*(z->m - y->m) >= (y->b - z->b)*(y->m - x->m);
    }
    void add(ll m, ll b) {
        auto y = insert({ m, b });
        y->succ = [=] { return next(y) == end() ? 0 : &*next(y); };
        if (bad(y)) { erase(y); return; }
        while (next(y) != end() && bad(next(y))) erase(next(y));
        while (y != begin() && bad(prev(y))) erase(prev(y));
    }
    ll query(ll x) {
        auto l = *lower_bound((Line) { x, is_query });
        return l.m * x + l.b;
    }
};
}

int test2() {
	LineContainer mh;
	const int K = 10;
	ll x[K], v[K];
	rep(it,0,100) {
		mh.clear();
		int N = rand() % 100000 + 1;
		rep(j,0,K) x[j] = rand() % 1000 - 500, v[j] = LLONG_MIN;
// cerr << "---" << endl;
// cerr << x << endl;
		rep(i,0,N) {
			ll k = rand() % 100000 - 50000;
			ll m = rand() % (1LL << 30) - (1LL << 29);
// cerr << k << ' ' << m << endl;
			mh.add((int)k, (int)m);
			rep(j,0,K) v[j] = max(v[j], k*x[j] + m);
		}
// cerr << mh.eval(x) << ' ' << v << endl;
// for(auto &li: mh) cerr << li.k << ' ' << li.m << ' ' << li.p << endl;
		rep(j,0,K)
			assert(mh.query(x[j]) == v[j]);
	}
	return 0;
}

// brute force with __int128, large coefficients, duplicates, sorted orders
void test3() {
	mt19937_64 rng(99);
	auto rnd = [&](ll lim) { return (ll)(rng() % (2 * (unsigned long long)lim + 1)) - lim; };
	rep(it,0,6000) {
		int N = (int)(rng() % 12) + 1;
		if (it % 100 == 0) N = 400;
		ll K = 1, M = 1, X = 1;
		switch (it % 6) {
			case 0: K = 3, M = 3, X = 5; break;
			case 1: K = 1000000000, M = (ll)1e18, X = 1000000000; break;
			case 2: K = 2, M = (ll)1e18, X = 1000000000; break;
			case 3: K = 1000000000, M = 2, X = 1000000000; break;
			case 4: K = (ll)1e6, M = (ll)1e12, X = (ll)1e6; break;
			case 5: K = (ll)1e18, M = (ll)1e18, X = 2; break;
		}
		vector<pair<ll, ll>> ls(N);
		for (auto& l : ls) l = {rnd(K), rnd(M)};
		int order = (int)(rng() % 4);
		if (order == 1) sort(all(ls));
		if (order == 2) sort(all(ls)), reverse(all(ls));
		if (order == 3) rep(i,1,N) if (rng() % 3 == 0) ls[i] = ls[rng() % i];
		LineContainer lc;
		rep(i,0,N) {
			lc.add(ls[i].first, ls[i].second);
			rep(q,0,N > 100 ? 2 : 6) {
				ll x = rnd(X);
				if (rng() % 8 == 0) x = rng() % 2 ? X : -X;
				__int128 best = (__int128)ls[0].first * x + ls[0].second;
				rep(j,1,i+1) best = max(best, (__int128)ls[j].first * x + ls[j].second);
				assert(best == lc.query(x));
			}
			// hull invariants: sorted slopes, strictly increasing breakpoints
			assert(prev(lc.end())->p == LineContainer::inf);
			for (auto a = lc.begin(), b = next(a); b != lc.end(); ++a, ++b)
				assert(a->k <= b->k && a->p < b->p);
		}
	}
}

volatile ll glob;
int ra() {
	static unsigned blah;
	blah *= 12311231;
	blah += 129481762;
	return blah >> 1;
}

int main() {
	LineContainer mh;
	other::HullDynamic mh2;
	rep(it,0,10000000) {
		assert(mh.empty() == mh2.empty());
		int r = ra() % 100;
		if (r < 10) mh.clear(), mh2.clear();
		else if (r < 50) {
			int k = ra() % 10 - 5;
			int m = ra() % 100 - 50;
			mh.add(k, m);
			mh2.add(k, m);
		}
		else if (!mh.empty()) {
			int x = ra() % 10 - 5;
			assert(mh.query(x) == mh2.query(x));
		}
	}
	test2();
	test3();
	cout<<"Tests passed!"<<endl;
}
