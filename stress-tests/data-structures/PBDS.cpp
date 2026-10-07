#include "../utilities/template.h"

#include "../../content/data-structures/PBDS.h"

// PBDS.h is a set of usage examples (order statistic tree, mergeable heap with
// modify, rope). Run the example itself and stress each container the way the
// example uses it against std:: oracles.

#ifdef CHECK_STD_NAMES
// Opt-in (fails to compile): after pasting PBDS.h, these std names are
// ambiguous because of `using namespace __gnu_pbds / __gnu_cxx`.
void stdNames() {
	priority_queue<int> q; // also __gnu_pbds::priority_queue
	q.push(1);
	vi v(3), w(3);
	copy_n(v.begin(), 2, w.begin()); // also __gnu_cxx::copy_n
}
#endif

mt19937_64 rng(2718);
int ri(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }

void testTree(int ops, int range) {
	Tree<int> t;
	vi v;
	rep(it,0,ops) {
		int x = ri(-range, range), op = ri(0, 5);
		auto pos = lower_bound(all(v), x);
		bool has = pos != v.end() && *pos == x;
		int idx = (int)(pos - v.begin());
		if (op <= 1) {
			assert(t.insert(x).second == !has);
			if (!has) v.insert(pos, x);
		} else if (op == 2) {
			assert(t.erase(x) == has);
			if (has) v.erase(pos);
		} else if (op == 3) assert((int)t.order_of_key(x) == idx);
		else if (op == 4 && sz(v)) {
			int k = ri(0, sz(v) - 1);
			assert(*t.find_by_order(k) == v[k]);
			assert(t.find_by_order(sz(v)) == t.end());
		} else { // join with a tree whose keys are all larger / all smaller
			Tree<int> o;
			int k = ri(0, 3), up = ri(0, 1);
			rep(i,0,k) {
				int y = up ? range + 1 + ri(0, 5) : -range - 1 - ri(0, 5);
				if (o.insert(y).second) v.push_back(y);
			}
			if (ri(0, 1)) t.join(o); else o.join(t), t.swap(o);
			sort(all(v));
			v.erase(unique(all(v)), v.end());
			// keep keys inside [-range, range] for the next join
			while (sz(v) && v.back() > range) t.erase(v.back()), v.pop_back();
			while (sz(v) && v[0] < -range) t.erase(v[0]), v.erase(v.begin());
		}
		assert(sz(t) == sz(v));
		if (it % 16 == 0) assert(equal(all(t), all(v)));
	}
}

void testHeap(int ops, int range) {
	const int K = 4;
	Heap<int> h[K];
	multiset<int> ms[K];
	// live handles: (heap index, iterator, current value)
	vector<tuple<int, Heap<int>::point_iterator, int>> hs;
	rep(it,0,ops) {
		int i = ri(0, K - 1), j = ri(0, K - 1), op = ri(0, 6), x = ri(-range, range);
		if (op <= 1) {
			hs.emplace_back(i, h[i].push(x), x);
			ms[i].insert(x);
		} else if (op == 2 && sz(ms[i])) {
			int top = h[i].top();
			assert(top == *ms[i].rbegin());
			h[i].pop();
			ms[i].erase(prev(ms[i].end()));
			// drop one handle with that value in heap i (equal values are interchangeable
			// for the oracle, but the popped iterator is unknown: drop them all)
			vector<tuple<int, Heap<int>::point_iterator, int>> keep;
			for (auto& t : hs) if (get<0>(t) != i || get<2>(t) != top) keep.push_back(t);
			hs = keep;
		} else if (op == 3 && i != j) { // meld j into i; iterators stay valid
			h[i].join(h[j]);
			assert(h[j].empty());
			ms[i].insert(all(ms[j])), ms[j].clear();
			for (auto& t : hs) if (get<0>(t) == j) get<0>(t) = i;
		} else if (op == 4 && sz(hs)) { // modify-key (both increase and decrease)
			int k = ri(0, sz(hs) - 1);
			auto& [hi, pit, val] = hs[k];
			assert(*pit == val);
			ms[hi].erase(ms[hi].find(val));
			h[hi].modify(pit, x);
			ms[hi].insert(val = x);
		} else if (op == 5 && sz(hs)) { // erase through the iterator
			int k = ri(0, sz(hs) - 1);
			auto [hi, pit, val] = hs[k];
			h[hi].erase(pit);
			ms[hi].erase(ms[hi].find(val));
			hs.erase(hs.begin() + k);
		}
		rep(k,0,K) {
			assert(h[k].size() == ms[k].size() && h[k].empty() == ms[k].empty());
			if (sz(ms[k])) assert(h[k].top() == *ms[k].rbegin());
		}
	}
	rep(k,0,K) while (sz(ms[k])) { // drain in sorted order
		assert(h[k].top() == *ms[k].rbegin());
		h[k].pop(), ms[k].erase(prev(ms[k].end()));
	}
}

void testRope(int ops, int maxLen) {
	rope<int> r;
	vi v;
	vector<rope<int>> oldR; vector<vi> oldV; // ropes are persistent by value
	rep(it,0,ops) {
		int op = ri(0, 6), n = sz(v);
		if (op == 0 && n < maxLen) {
			int x = ri(-9, 9);
			r.push_back(x), v.push_back(x);
		} else if (op == 1 && n) {
			int i = ri(0, n - 1), x = ri(-9, 9);
			r.mutable_reference_at(i) = x, v[i] = x;
		} else if (op == 2) { // erase(pos, len)
			int l = ri(0, n), len = ri(0, n - l);
			r.erase(l, len), v.erase(v.begin() + l, v.begin() + l + len);
		} else if (op == 3 && n) { // move substr(l, len) to position p
			int l = ri(0, n - 1), len = ri(1, n - l);
			rope<int> cur = r.substr(l, len);
			vi cv(v.begin() + l, v.begin() + l + len);
			r.erase(l, len), v.erase(v.begin() + l, v.begin() + l + len);
			int p = ri(0, sz(v));
			if (ri(0, 1)) r.insert(r.mutable_begin() + p, cur); else r.insert(p, cur);
			v.insert(v.begin() + p, all(cv));
		} else if (op == 4 && n && 2 * n <= maxLen) { // insert a copy of a substring
			int l = ri(0, n - 1), len = ri(1, n - l), p = ri(0, n);
			rope<int> cur = r.substr(l, len);
			vi cv(v.begin() + l, v.begin() + l + len);
			r.insert(p, cur), v.insert(v.begin() + p, all(cv));
		} else if (op == 5) {
			if (sz(oldR) < 8) oldR.push_back(r), oldV.push_back(v);
			else { int k = ri(0, 7); r = oldR[k], v = oldV[k]; }
		} else if (n) {
			int i = ri(0, n - 1);
			assert(r[i] == v[i] && r.at(i) == v[i]);
		}
		assert(sz(r) == sz(v));
		if (it % 8 == 0) {
			assert(equal(r.begin(), r.end(), v.begin()));
			rep(k,0,sz(oldR)) assert(sz(oldR[k]) == sz(oldV[k]) &&
				equal(oldR[k].begin(), oldR[k].end(), oldV[k].begin()));
		}
	}
}

// The rope part of pbds() in the header, step by step, with the states its
// comments claim.
void ropeExample() {
	auto is = [](const rope<int>& r, vi e) {
		return sz(r) == sz(e) && equal(r.begin(), r.end(), e.begin());
	};
	int n = 3;
	rope<int> v(n, 0);
	for (int i=0; i<n; i++) v.mutable_reference_at(i) = i + 1;
	assert(is(v, {1, 2, 3}));
	for (int i=0; i<n; i++) v.push_back(i + n + 1);
	assert(is(v, {1, 2, 3, 4, 5, 6}));
	int l=1, r=3;
	rope<int> cur = v.substr(l, r-l+1);
	assert(is(cur, {2, 3, 4}));
	v.erase(l, r-l+1);
	assert(is(v, {1, 5, 6}));
	rope<int> a = v, b = v, c = v;
	a.insert(a.mutable_begin(), cur);
	assert(is(a, {2, 3, 4, 1, 5, 6})); // "to start"
	b.insert(b.mutable_begin() + 2, cur);
	assert(is(b, {1, 5, 2, 3, 4, 6})); // "to TWO AFTER start"
	c.insert(c.mutable_begin() + 1, cur);
	assert(is(c, {1, 2, 3, 4, 5, 6})); // ONE AFTER start
	// `v.insert(v.mutable_reference_at(0), cur)` converts the element v[0] to a
	// position: it inserts at index v[0], which is 1 only because v[0] == 1.
	rope<int> d = v, e = v;
	d.insert(d.mutable_reference_at(0), cur);
	assert(is(d, {1, 2, 3, 4, 5, 6}));
	e.mutable_reference_at(0) = 3;
	e.insert(e.mutable_reference_at(0), cur);
	assert(is(e, {3, 5, 6, 2, 3, 4}));
}

int main(int argc, char**) {
	if (argc > 1) { // benchmark mode: ./a.out bench
		auto now = [] { return chrono::steady_clock::now(); };
		auto secs = [&](auto t0) { return chrono::duration<double>(now() - t0).count(); };
		const int N = 1000000;
		{
			Tree<int> t;
			auto t0 = now();
			ll s = 0;
			rep(i,0,N) t.insert((int)rng());
			rep(i,0,N) s += t.order_of_key((int)rng()) + *t.find_by_order(rng() % sz(t));
			cerr << "Tree: 1e6 insert + 1e6 order_of_key + 1e6 find_by_order: " << secs(t0) << " s (" << s << ")\n";
		}
		{
			Heap<int> h;
			std::priority_queue<int> q;
			vector<Heap<int>::point_iterator> its(N);
			auto t0 = now();
			rep(i,0,N) its[i] = h.push((int)rng());
			rep(i,0,N) h.modify(its[rng() % N], (int)rng());
			ll s = 0;
			rep(i,0,N) s += h.top(), h.pop();
			double th = secs(t0); t0 = now();
			rep(i,0,N) q.push((int)rng());
			rep(i,0,N) s += q.top(), q.pop();
			cerr << "Heap: 1e6 push + 1e6 modify + 1e6 pop: " << th << " s; std::priority_queue 1e6 push+pop: "
				<< secs(t0) << " s (" << s << ")\n";
			vector<Heap<int>> hs(N); // small-to-large free melding
			rep(i,0,N) hs[i].push(i);
			t0 = now();
			rep(i,1,N) hs[0].join(hs[i]);
			cerr << "Heap: 1e6 joins: " << secs(t0) << " s\n";
			while (!hs[0].empty()) hs[0].pop(); // a deep heap must not be destroyed as is
		}
		for (int n : {100000, 1000000}) {
			rope<int> r;
			auto t0 = now();
			rep(i,0,n) r.push_back(i);
			double tp = secs(t0); t0 = now();
			int Q = 100000;
			rep(i,0,Q) { // cut [l, l+len) and paste at the front
				int l = ri(0, n - 1), len = ri(1, n - l);
				rope<int> cur = r.substr(l, len);
				r.erase(l, len);
				r.insert(r.mutable_begin(), cur);
			}
			double tc = secs(t0); t0 = now();
			ll s = 0;
			rep(i,0,Q) s += r[ri(0, n - 1)];
			double ta = secs(t0); t0 = now();
			rep(i,0,Q) r.mutable_reference_at(ri(0, n - 1)) = i;
			cerr << "rope n=" << n << ": push_back x n " << tp << " s, 1e5 cut+paste " << tc
				<< " s, 1e5 random reads " << ta << " s, 1e5 mutable_reference_at " << secs(t0)
				<< " s (" << s << ")\n";
		}
		return 0;
	}
	pbds();
	ropeExample();
	rep(it,0,20000) testTree(60, 4);
	rep(it,0,1000) testTree(2000, 60);
	testTree(20000, 1000000000);
	rep(it,0,20000) testHeap(60, 3);
	rep(it,0,1000) testHeap(2000, 1000000000);
	rep(it,0,20000) testRope(60, 8);
	rep(it,0,1000) testRope(2000, 200);
	testRope(200000, 5000);
#ifdef HEAP_DEEP_DESTRUCT
	// Opt-in (segfaults with the default 8 MB stack; `make test` raises the
	// limit): pushing increasing keys builds a pairing heap that is a path, and
	// its destructor / copy constructor recurse once per element.
	{
		Heap<int> h;
		rep(i,0,1000000) h.push(i);
	}
#endif
	cout << "Tests passed!" << endl;
}
