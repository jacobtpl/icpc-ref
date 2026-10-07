#include "../utilities/template.h"

#include "../../content/data-structures/OrderStatisticTree.h"

mt19937_64 rng(31337);
ll rl(ll a, ll b) { return uniform_int_distribution<ll>(a, b)(rng); }

template<class T, class G> void test(int ops, G gen) {
	Tree<T> t;
	vector<T> v; // sorted oracle
	rep(it,0,ops) {
		T x = gen();
		auto pos = lower_bound(all(v), x);
		bool has = pos != v.end() && *pos == x;
		int idx = (int)(pos - v.begin()), op = (int)rl(0, 6);
		if (op <= 1) {
			auto r = t.insert(x);
			assert(r.second == !has && *r.first == x); // a set: duplicates rejected
			if (!has) v.insert(pos, x);
		} else if (op == 2) {
			assert(t.erase(x) == has);
			if (has) v.erase(pos);
		} else if (op == 3) {
			assert((int)t.order_of_key(x) == idx);
		} else if (op == 4) {
			int k = (int)rl(0, sz(v) + 1);
			auto f = t.find_by_order(k);
			if (k < sz(v)) assert(f != t.end() && *f == v[k]);
			else assert(f == t.end());
		} else if (op == 5) {
			auto lb = t.lower_bound(x), ub = t.upper_bound(x);
			assert(lb == t.end() ? idx == sz(v) : *lb == v[idx]);
			int j = idx + has;
			assert(ub == t.end() ? j == sz(v) : *ub == v[j]);
		} else if (sz(v)) { // split at a key, check both parts, join back
			Tree<T> hi;
			t.split(x, hi); // t keeps keys <= x
			int j = idx + has;
			assert(sz(t) == j && sz(hi) == sz(v) - j);
			if (j < sz(v)) assert(*hi.begin() == v[j] && (int)hi.order_of_key(v.back()) == sz(v) - j - 1);
			if (rl(0, 1)) t.join(hi); else hi.join(t), t.swap(hi);
			assert(hi.empty());
		}
		assert(sz(t) == sz(v));
		if (it % 16 == 0) assert(equal(all(t), all(v)));
	}
}

int main(int argc, char**) {
	if (argc > 1) { // benchmark mode: ./a.out bench
		auto now = [] { return chrono::steady_clock::now(); };
		auto secs = [&](auto t0) { return chrono::duration<double>(now() - t0).count(); };
		const int N = 1000000;
		for (int sorted : {0, 1}) {
			vi xs(N);
			rep(i,0,N) xs[i] = sorted ? i : (int)rng();
			Tree<int> t;
			set<int> s;
			auto t0 = now();
			for (int x : xs) t.insert(x);
			double ti = secs(t0); t0 = now();
			ll h = 0;
			for (int x : xs) h += t.order_of_key(x);
			double to = secs(t0); t0 = now();
			rep(i,0,N) h += *t.find_by_order(rng() % sz(t));
			double tf = secs(t0); t0 = now();
			for (int x : xs) t.erase(x);
			double te = secs(t0); t0 = now();
			for (int x : xs) s.insert(x);
			double ts = secs(t0);
			cerr << (sorted ? "sorted" : "random") << " keys, N=1e6: insert " << ti
				<< " s, order_of_key " << to << " s, find_by_order " << tf << " s, erase " << te
				<< " s; std::set insert " << ts << " s (" << h << ")\n";
		}
		return 0;
	}
	example();
	rep(it,0,3000) test<int>(60, [] { return (int)rl(-4, 4); });
	rep(it,0,300) test<int>(1000, [] { return (int)rl(-40, 40); });
	rep(it,0,300) test<int>(1000, [] { return rl(0, 3) ? (int)rl(INT_MIN, INT_MAX) : rl(0, 1) ? INT_MIN : INT_MAX; });
	rep(it,0,300) test<pair<ll, int>>(1000, [] { return make_pair(rl(-2, 2) * (LLONG_MAX / 2), (int)rl(0, 3)); });
	test<int>(300000, [] { return (int)rl(0, 5000); });
	{ // empty tree, and the documented join precondition (disjoint key ranges)
		Tree<int> a, b;
		assert(a.order_of_key(5) == 0 && a.find_by_order(0) == a.end());
		a.join(b);
		assert(a.empty());
		a.insert(1), a.insert(5), b.insert(3);
		bool thrown = 0;
		try { a.join(b); } catch (join_error&) { thrown = 1; }
		assert(thrown && sz(a) == 2 && sz(b) == 1);
	}
	cout << "Tests passed!" << endl;
}
