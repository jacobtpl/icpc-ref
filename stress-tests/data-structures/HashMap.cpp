#include "../utilities/template.h"

#include "../../content/data-structures/HashMap.h"

typedef __gnu_pbds::gp_hash_table<ll, int, chash> HT;

template<class H>
void check(H& t, map<ll, int>& m) {
	assert(t.size() == m.size() && t.empty() == m.empty());
	map<ll, int> seen;
	for (auto& p : t) assert(seen.insert(p).second);
	assert(seen == m);
}

template<class H, class F>
void run(H& t, int ops, int checkEvery, F key, mt19937_64& rng) {
	map<ll, int> m;
	t.clear();
	rep(it,0,ops) {
		ll k = key();
		int v = (int)rng(), op = (int)(rng() % 100);
		if (op < 30) t[k] = v, m[k] = v;
		else if (op < 40) t[k] += v % 1000, m[k] += v % 1000;
		else if (op < 50) {
			auto a = t.insert({k, v}); auto b = m.insert({k, v});
			assert(a.second == b.second && a.first->second == b.first->second);
		}
		else if (op < 70) assert(t.erase(k) == (m.erase(k) > 0));
		else if (op < 99) {
			auto a = t.find(k); auto b = m.find(k);
			assert((a == t.end()) == (b == m.end()));
			if (b != m.end()) assert(a->first == k && a->second == b->second);
		}
		else if (rng() % 20 == 0) t.clear(), m.clear();
		if (it % checkEvery == 0 || it == ops - 1) check(t, m);
	}
}

int main() {
	mt19937_64 rng(2024);
	assert(h.empty() && h.find(0) == h.end());
	const ll edge[] = {0, 1, -1, LLONG_MIN, LLONG_MAX, LLONG_MIN + 1, LLONG_MAX - 1,
		1LL << 16, 1LL << 32, 1LL << 48, -(1LL << 32), INT_MIN, INT_MAX};
	uint64_t inv = 1, C = chash().C; // inverse of C mod 2^64
	rep(i,0,6) inv *= 2 - C * inv;
	assert(inv * C == 1);
	vector<function<ll()>> gens = {
		[&]() { return (ll)(rng() % 20) - 10; },
		[&]() { return (ll)(rng() % 3000); },
		[&]() { return (ll)rng(); },
		[&]() { return (ll)((uint64_t)edge[rng() % 13] + (rng() % 3 == 0)); },
		[&]() { return (ll)(rng() % 300) << (16 * (rng() % 4)); }, // equal low bits
		[&]() { return (ll)((rng() % 300) * inv); }, // all hash to the same bucket
	};
	rep(g,0,sz(gens)) {
		rep(it,0,300) run(h, 300, 60, gens[g], rng); // iterating h is slow: 1 << 16 buckets
		rep(it,0,3) run(h, 100000, 5000, gens[g], rng);
		// tables that start small and have to grow
		rep(it,0,100) { HT t; run(t, 500, 1, gens[g], rng); }
		rep(it,0,100) { HT t({},{},{},{},{(size_t)1 << (rng() % 5)}); run(t, 500, 1, gens[g], rng); }
		rep(it,0,2) { HT t({},{},{},{},{2}); run(t, 100000, 5000, gens[g], rng); }
	}
	{ // growth far beyond the initial 1 << 16 buckets
		h.clear();
		unordered_map<ll, int> m;
		rep(i,0,1000000) {
			ll k = i % 2 ? (ll)rng() : i;
			h[k] += i, m[k] += i;
		}
		assert(h.size() == m.size());
		for (auto& p : m) assert(h.find(p.first)->second == p.second);
		for (auto& p : m) if (p.first & 1) assert(h.erase(p.first));
		for (auto& p : m) assert((h.find(p.first) == h.end()) == (p.first & 1));
	}
	cout<<"Tests passed!"<<endl;
}
