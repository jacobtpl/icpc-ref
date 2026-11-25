#include "../utilities/template.h"

#include "../../content/data-structures/UnionFindPotential.h"

namespace X {
#define PotentialUF XorUF
struct PotentialUF {
	typedef ll T;
	static T add(T a, T b) { return a ^ b; }
	static T neg(T a) { return a; }

	vi e; vector<T> w; vector<pii> st;
	PotentialUF(int n) : e(n, -1), w(n) {}
	pair<int, T> find(int x) {
		T d{};
		for (; e[x] >= 0; x = e[x]) d = add(d, w[x]);
		return {x, d};
	}
	bool same(int a, int b) { return find(a).first == find(b).first; }
	T diff(int a, int b) {
		return add(find(a).second, neg(find(b).second));
	}
	int join(int a, int b, T v) {
		auto [ra, da] = find(a); auto [rb, db] = find(b);
		if (ra == rb) return add(da, neg(db)) == v ? 0 : -1;
		T x = add(add(v, db), neg(da));
		if (e[ra] < e[rb]) swap(ra, rb), x = neg(x);
		st.push_back({ra, e[ra]});
		e[rb] += e[ra]; e[ra] = rb; w[ra] = x;
		return 1;
	}
	int time() { return sz(st); }
	void rollback(int t) {
		for (; sz(st) > t; st.pop_back()) {
			auto [c, old] = st.back();
			e[e[c]] -= old; e[c] = old;
		}
	}
};
#undef PotentialUF
}

// Brute force: constraints as edges; BFS for the implied difference.
template<class U, class Add, class Neg>
void test(Add add, Neg neg, int iters) {
	rep(it,0,iters) {
		int n = rand() % 8 + 1;
		U uf(n);
		vector<tuple<int,int,ll>> cons; // accepted constraints
		vector<pair<int,int>> marks; // (uf.time(), sz(cons))
		auto brute = [&](int a, int b) -> pair<bool, ll> {
			vector<ll> d(n); vi seen(n); seen[a] = 1; // d[x] = p(a) - p(x)
			queue<int> q; q.push(a);
			while (!q.empty()) {
				int x = q.front(); q.pop();
				for (auto [u, v, w] : cons) { // p(u) - p(v) = w
					if (u == x && !seen[v]) seen[v] = 1, d[v] = add(d[x], w), q.push(v);
					if (v == x && !seen[u]) seen[u] = 1, d[u] = add(d[x], neg(w)), q.push(u);
				}
			}
			return {seen[b], d[b]};
		};
		rep(op,0,40) {
			int a = rand() % n, b = rand() % n, t = rand() % 10;
			if (t < 5) {
				ll v = rand() % 7 - 3;
				if (rand() % 2) { auto [ok, d] = brute(a, b); if (ok) v = d; }
				if (a == b && rand() % 2) v = 0;
				auto [ok, d] = brute(a, b);
				int exp = !ok ? 1 : d == v ? 0 : -1;
				// stored constraint is p(a) - p(b) = v; brute d = p(a) - p(b)
				int got = uf.join(a, b, v);
				assert(got == exp);
				if (got == 1) cons.emplace_back(a, b, v);
			} else if (t < 8) {
				auto [ok, d] = brute(a, b);
				assert(uf.same(a, b) == ok);
				if (ok) assert(uf.diff(a, b) == d);
				int cnt = 0;
				rep(x,0,n) cnt += brute(a, x).first;
				assert(uf.find(a).first == uf.find(uf.find(a).first).first);
				(void)cnt;
			} else if (t == 8) {
				marks.push_back({uf.time(), sz(cons)});
			} else if (!marks.empty()) {
				int k = rand() % sz(marks);
				uf.rollback(marks[k].first);
				cons.resize(marks[k].second);
				marks.resize(k);
			}
		}
	}
}

int main() {
	test<PotentialUF>([](ll a, ll b) { return a + b; }, [](ll a) { return -a; }, 20000);
	test<X::XorUF>([](ll a, ll b) { return a ^ b; }, [](ll a) { return a; }, 20000);
	// sizes
	rep(it,0,2000) {
		int n = rand() % 10 + 1; PotentialUF uf(n); vi c(n); iota(all(c), 0);
		rep(op,0,15) {
			int a = rand() % n, b = rand() % n;
			if (uf.join(a, b, rand() % 5) == 1) { int x = c[a], y = c[b]; for (int &z : c) if (z == x) z = y; }
			rep(i,0,n) assert(uf.size(i) == (int)count(all(c), c[i]));
		}
		uf.rollback(0);
		rep(i,0,n) assert(uf.size(i) == 1);
	}
	cout<<"Tests passed!"<<endl;
}
