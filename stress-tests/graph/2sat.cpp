#include "../utilities/template.h"

#include "../../content/graph/2sat.h"

int main1() {
	const int N = 100000, M = 10000000;
	// Random constraints, unsolvable
	{
		TwoSat ts(N);
		rep(i,0,M) {
			int r = rand();
			int s = r;
			r >>= 2;
			int a = r % N;
			r >>= 5;
			int b = r % N;
			if (a == b) continue;
			ts.either(a ^ (s&1 ? 0 : -1), b ^ (s&2 ? 0 : -1));
		}
		assert(ts.solve() == 0);
	}
	// Random solvable instance
	{
		vector<bool> v(N);
		rep(i,0,N) v[i] = rand() & (1 << 20);
		TwoSat ts(N);
		rep(i,0,M) {
			int r = rand();
			int s = r;
			r >>= 2;
			int a = r % N;
			r >>= 5;
			int b = r % N;
			if (a == b) continue;
			ts.either(a ^ (v[a] ? 0 : -1), b ^ (s&1 ? 0 : -1));
		}
		assert(ts.solve() == 1);
	}
	return 0;
}

int main2() {
	int N = 4;
	TwoSat ts(N);
	ts.either(0,1);
	ts.either(0,~1);
	ts.either(~2,~3);
	assert(ts.solve()==1);
	assert(ts.values == vi({1, 1, 0, 0}));
	return 0;
}

// Exhaustive check against brute force on tiny instances, including
// either/setValue/atMostOne with repeated and contradictory literals.
int ra();
void bruteTest() {
	rep(it,0,200000) {
		int N = ra() % 6 + 1, M = ra() % 12;
		TwoSat ts(N);
		vector<pii> cl;
		vector<vi> atm;
		auto lit = [&]() { int a = ra() % N; return ra() % 2 ? a : ~a; };
		rep(i,0,M) {
			int t = ra() % 8;
			if (t < 5) {
				int a = lit(), b = lit();
				ts.either(a, b);
				cl.emplace_back(a, b);
			} else if (t < 6) {
				int a = lit();
				ts.setValue(a);
				cl.emplace_back(a, a);
			} else {
				vi li(ra() % 5);
				for (int& x : li) x = lit();
				ts.atMostOne(li);
				atm.push_back(li);
			}
		}
		auto ok = [&](int mask) {
			auto val = [&](int x) { return x >= 0 ? mask >> x & 1 : !(mask >> ~x & 1); };
			for (auto [a, b] : cl) if (!val(a) && !val(b)) return false;
			for (auto& li : atm) {
				int c = 0;
				for (int x : li) c += val(x);
				if (c > 1) return false;
			}
			return true;
		};
		bool sat = 0;
		rep(mask,0,1<<N) if (ok(mask)) sat = 1;
		bool got = ts.solve();
		assert(got == sat);
		if (got) {
			int mask = 0;
			rep(i,0,N) {
				assert(ts.values[i] == 0 || ts.values[i] == 1);
				mask |= ts.values[i] << i;
			}
			assert(ok(mask));
			assert(ts.solve()); // solving twice is fine
		}
	}
	TwoSat empty;
	assert(empty.solve() && empty.values.empty());
	int a = empty.addVar();
	empty.setValue(~a);
	assert(empty.solve() && empty.values == vi({0}));
	empty.setValue(a);
	assert(!empty.solve());
}

int ra() {
	static unsigned X;
	X *= 1283611;
	X += 123;
	return X >> 1;
}

// Test at_most_one
int main() {
	main1();
	main2();
	bruteTest();
	const int N = 100, M = 400;
	rep(it,0,100) {
		vector<bool> v(N);
		rep(i,0,N) v[i] = ra() & (1 << 20);
		TwoSat ts(N);
		vector<vi> atm;
		vi r;
		rep(i,0,M) {
			if (ra()%100 < 5) {
				int r = ra();
				int s = r;
				r >>= 2;
				int a = r % N;
				r >>= 5;
				int b = r % N;
				if (a == b) continue;
				ts.either(v[a] ? a : ~a, (s&1) ? b : ~b);
			} else {
				int k = ra() % 4 + 1;
				r.clear();
				rep(ki,0,k-1) {
					int a = ra() % N;
					r.push_back(v[a] ? ~a : a);
				}
				r.push_back(ra() % (2*N) - N);
				random_shuffle(all(r), [](int x) { return ra() % x; });
				ts.atMostOne(r);
				atm.push_back(r);
			}
		}
		assert(ts.solve());
		int to = 0;
		rep(i,0,N) to += (ts.values[i] == v[i]);
		for(auto &r: atm) {
			int co = 0;
			for(auto &x: r) co += (ts.values[max(x, ~x)] == (x >= 0));
			assert(co <= 1);
		}
	}
	cout<<"Tests passed!"<<endl;
	return 0;
}
